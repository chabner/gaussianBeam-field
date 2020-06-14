#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238

__constant__ double box_min[3];
__constant__ double box_max[3];
__constant__ double sigt; // sigt/2

__constant__ int32_t dirDimNum; // which dim the directions are? 2/3/4

__constant__ int32_t uDimProd[5];
__constant__ int32_t lDim[5];
__constant__ int32_t vDim[5];
__constant__ int32_t mixturesNum;

__constant__ double lMixtureAlpha[32];


__forceinline __device__ double2 operator +(double2 a, double2 b)
{
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    return a;
}

__forceinline __device__ double2 operator +(double2 a, double b)
{
    a.x = a.x + b;
    return a;
}

__forceinline __device__ double2 operator /(double2 a, double2 b)
{
    double2 c;
    double denominator = 1 / (fma(b.x, b.x, b.y * b.y));
    c.x = (fma(a.x, b.x, a.y * b.y)) * denominator;
    c.y = (fma(a.y, b.x, -a.x * b.y)) * denominator;
    return c;
}

__forceinline __device__ double2 operator *(double a, double2 b)
{
    b.x = b.x * a;
    b.y = b.y * a;
    return b;
}

__forceinline __device__ double2 operator *(double2 a, double2 b)
{
    register double2 tmpVar = a;

    a.x = fma(a.x, b.x, -a.y * b.y);
    a.y = fma(tmpVar.x, b.y, tmpVar.y * b.x);
    return a;
}

__forceinline __device__ double2 complexSquare(double2 a)
{
    register double tmpVar = a.x;
    a.x = (a.x - a.y) * (a.x + a.y);
    a.y = 2 * tmpVar * a.y;
    return a;
}

__forceinline __device__ double2 complexConj(double2 a)
{
    a.y = -1 * a.y;
    return a;
}

__device__ double2 complexSqrt(double2 a)
{
    double r, absrz;
    r = hypot(a.x, a.y);
    absrz = sqrt(r) * rhypot(r + a.x, a.y);
    // absrz = sqrt(r) / (fma(sqrt((r + a.x)),(r + a.x),a.y * a.y));
    a.x = absrz * (a.x + r);
    a.y = absrz * a.y;
    return a;
}

__forceinline __device__ double2 complexExponent(double2 a)
{
    register double expr, sina, cosa;
    expr = exp(a.x);
    sincos(a.y, &sina, &cosa);
    a.x = expr * cosa;
    a.y = expr * sina;
    return a;
}

__forceinline __device__ void complexIncrease(double2* address, double2 val)
{
    //(*address).x += val.x;
    //(*address).y += val.y;

    atomicAdd(&(address->x), val.x);
    atomicAdd(&(address->y), val.y);
}

__global__ void MSintegrateMult_correlation(double2* u, const double2* el, const double* dirv,
    const double* x, const double2* randomPhase, const bool* activatedPaths,
    const double2* vApertureMu1, const double2* vApertureMu2, const double2* vApertureMu3, const double2* vApertureC,
    const double* mixtureMu1, const double* mixtureMu2, const double* mixtureMu3, const double* mixtureC)
{
    int bIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (bIdx < uDimProd[4])
    {
        int uDimIdx[5];
        uDimIdx[4] = bIdx / (uDimProd[3]);
        uDimIdx[3] = 0;
        uDimIdx[2] = (bIdx - uDimIdx[4] * (uDimProd[3]) - uDimIdx[3] * (uDimProd[2])) / uDimProd[1];
        uDimIdx[1] = (bIdx - uDimIdx[4] * (uDimProd[3]) - uDimIdx[3] * (uDimProd[2]) - uDimIdx[2] * (uDimProd[1])) / uDimProd[0];
        uDimIdx[0] = bIdx % uDimProd[0];

        if (activatedPaths[uDimIdx[4]])
        {
            int uIdx = uDimIdx[0] + uDimIdx[1] * uDimProd[0] + uDimIdx[2] * uDimProd[1];

            int lIdx = (lDim[1] > 1)* uDimIdx[0] + lDim[1] * (lDim[2] > 1)* uDimIdx[1] + lDim[1] * lDim[2] * (lDim[3] > 1)* uDimIdx[2] +
                lDim[1] * lDim[2] * lDim[3] * uDimIdx[3] + lDim[1] * lDim[2] * lDim[3] * lDim[4] * uDimIdx[4];

            int vIdx = (vDim[1] > 1)* uDimIdx[0] + vDim[1] * (vDim[2] > 1)* uDimIdx[1] + vDim[1] * vDim[2] * (vDim[3] > 1)* uDimIdx[2] +
                vDim[1] * vDim[2] * vDim[3] * uDimIdx[3];

            int dirIdx = uDimIdx[dirDimNum - 2] + vDim[dirDimNum - 1] * uDimIdx[3];

            double2 sqrtMu_1, sqrtMu_2, randPhase, u1, u2;
            double2 vThroughputMu_1[3], vThroughputC_1, vThroughputMu_2[3], vThroughputC_2;

            double bd1, bd2, bd3, dz, alpha;

            const double* currX = x + uDimIdx[4] * 3;
            int currMixtureIdx_x;

            /* Aperture 1 */
            double2 elVal_1 = el[lIdx];
            const double* v = dirv + dirIdx * 3;

            // Cube dist
            if (v[0] <= 0)
            {
                bd1 = (currX[0] - box_min[0]) / abs(v[0]);
            }
            else
            {
                bd1 = (box_max[0] - currX[0]) / abs(v[0]);
            }

            if (v[1] <= 0)
            {
                bd2 = (currX[1] - box_min[1]) / abs(v[1]);
            }
            else
            {
                bd2 = (box_max[1] - currX[1]) / abs(v[1]);
            }

            if (v[2] <= 0)
            {
                bd3 = (currX[2] - box_min[2]) / abs(v[2]);
            }
            else
            {
                bd3 = (box_max[2] - currX[2]) / abs(v[2]);
            }

            dz = abs(fmin(fmin(bd1, bd2), bd3));

            // Throughput
            vThroughputMu_1[0] = vApertureMu1[vIdx];
            vThroughputMu_1[0].y -= 2 * PI * currX[0];

            vThroughputMu_1[1] = vApertureMu2[vIdx];
            vThroughputMu_1[1].y -= 2 * PI * currX[1];

            vThroughputMu_1[2] = vApertureMu3[vIdx];
            vThroughputMu_1[2].y -= 2 * PI * currX[2];

            vThroughputC_1 = vApertureC[vIdx] + (-sigt * dz);

            /* Aperture 2 */
            uDimIdx[3] = 1;

            lIdx = lIdx + lDim[0] * lDim[1] * lDim[2] * lDim[3];

            vIdx = vIdx + vDim[1] * vDim[2] * vDim[3];

            dirIdx = uDimIdx[dirDimNum - 2] + vDim[dirDimNum - 1];

            double2 elVal_2 = el[lIdx];
            v = dirv + dirIdx * 3;

            // Cube dist
            if (v[0] <= 0)
            {
                bd1 = (currX[0] - box_min[0]) / abs(v[0]);
            }
            else
            {
                bd1 = (box_max[0] - currX[0]) / abs(v[0]);
            }

            if (v[1] <= 0)
            {
                bd2 = (currX[1] - box_min[1]) / abs(v[1]);
            }
            else
            {
                bd2 = (box_max[1] - currX[1]) / abs(v[1]);
            }

            if (v[2] <= 0)
            {
                bd3 = (currX[2] - box_min[2]) / abs(v[2]);
            }
            else
            {
                bd3 = (box_max[2] - currX[2]) / abs(v[2]);
            }

            dz = abs(fmin(fmin(bd1, bd2), bd3));

            // Throughput
            vThroughputMu_2[0] = vApertureMu1[vIdx];
            vThroughputMu_2[0].y -= 2 * PI * currX[0];

            vThroughputMu_2[1] = vApertureMu2[vIdx];
            vThroughputMu_2[1].y -= 2 * PI * currX[1];

            vThroughputMu_2[2] = vApertureMu3[vIdx];
            vThroughputMu_2[2].y -= 2 * PI * currX[2];

            vThroughputC_2 = vApertureC[vIdx] + (-sigt * dz);

            // Integrate Mult
            currMixtureIdx_x = mixturesNum * uDimIdx[4];
            randPhase = randomPhase[uDimIdx[4]];

            u1.x = 0;
            u1.y = 0;
            u2 = u1;

            for (int mixtureIdx = 0; mixtureIdx < mixturesNum; mixtureIdx++)
            {
                sqrtMu_1 = complexSqrt(complexSquare(vThroughputMu_1[0] + mixtureMu1[currMixtureIdx_x]) +
                    complexSquare(vThroughputMu_1[1] + mixtureMu2[currMixtureIdx_x]) +
                    complexSquare(vThroughputMu_1[2] + mixtureMu3[currMixtureIdx_x]));

                sqrtMu_2 = complexSqrt(complexSquare(vThroughputMu_2[0] + mixtureMu1[currMixtureIdx_x]) +
                    complexSquare(vThroughputMu_2[1] + mixtureMu2[currMixtureIdx_x]) +
                    complexSquare(vThroughputMu_2[2] + mixtureMu3[currMixtureIdx_x]));

                
                alpha = lMixtureAlpha[mixtureIdx];

                u1 = u1 + alpha * randPhase * elVal_1 * complexExponent(vThroughputC_1 + mixtureC[mixtureIdx] + sqrtMu_1) / sqrtMu_1;
                u2 = u2 + alpha * randPhase * elVal_2 * complexExponent(vThroughputC_2 + mixtureC[mixtureIdx] + sqrtMu_2) / sqrtMu_2;

                currMixtureIdx_x++;
            }

            complexIncrease(u + uIdx, u1 * complexConj(u2));
        }
    }
}

int main() {
    return 0;

}