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
    (*address).x += val.x;
    (*address).y += val.y;
}

__global__ void integrateMult(double2* u, const double* dirv, const double* x, const double2* randomPhase,
    const double2* lMixtureMu1, const double2* lMixtureMu2, const double2* lMixtureMu3, const double2* lMixtureC,
    const double2* vApertureMu1, const double2* vApertureMu2, const double2* vApertureMu3, const double2* vApertureC)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < uDimProd[4])
    {
        int uDimIdx[5];
        uDimIdx[4] = uIdx / (uDimProd[3]);
        uDimIdx[3] = (uIdx - uDimIdx[4] * (uDimProd[3])) / uDimProd[2];
        uDimIdx[2] = (uIdx - uDimIdx[4] * (uDimProd[3]) - uDimIdx[3] * (uDimProd[2])) / uDimProd[1];
        uDimIdx[1] = (uIdx - uDimIdx[4] * (uDimProd[3]) - uDimIdx[3] * (uDimProd[2]) - uDimIdx[2] * (uDimProd[1])) / uDimProd[0];
        uDimIdx[0] = uIdx % uDimProd[0];

        int lIdx = lDim[0] * (lDim[1] > 1)* uDimIdx[0] + lDim[0] * lDim[1] * (lDim[2] > 1)* uDimIdx[1] + lDim[0] * lDim[1] * lDim[2] * (lDim[3] > 1)* uDimIdx[2] +
            lDim[0] * lDim[1] * lDim[2] * lDim[3] * uDimIdx[3] + lDim[0] * lDim[1] * lDim[2] * lDim[3] * lDim[4] * uDimIdx[4];

        int vIdx = (vDim[1] > 1)* uDimIdx[0] + vDim[1] * (vDim[2] > 1)* uDimIdx[1] + vDim[1] * vDim[2] * (vDim[3] > 1)* uDimIdx[2] +
            vDim[1] * vDim[2] * vDim[3] * uDimIdx[3];

        int dirIdx = uDimIdx[dirDimNum - 2] + vDim[dirDimNum - 1] * uDimIdx[3];

        double2 sqrtMu;
        int currlIdx;
        const double* currX = x + uDimIdx[4] * 3;
        double2 vThroughputMu[3], vThroughputC;

        const double* v = dirv + dirIdx * 3;

        double bd1, bd2, bd3, dz;

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

        dz = fmin(fmin(bd1, bd2), bd3);

        // Throughput
        vThroughputMu[0] = vApertureMu1[vIdx];
        vThroughputMu[0].y -= 2 * PI * currX[0];

        vThroughputMu[1] = vApertureMu2[vIdx];
        vThroughputMu[1].y -= 2 * PI * currX[1];

        vThroughputMu[2] = vApertureMu3[vIdx];
        vThroughputMu[2].y -= 2 * PI * currX[2];

        vThroughputC = vApertureC[vIdx] + (-sigt * dz);

        // Integrate Mult
        for (int mixtureIdx = 0; mixtureIdx < lDim[0]; mixtureIdx++)
        {
            currlIdx = lIdx + mixtureIdx;
            sqrtMu = complexSqrt(complexSquare(lMixtureMu1[currlIdx] + vThroughputMu[0]) +
                complexSquare(lMixtureMu2[currlIdx] + vThroughputMu[1]) +
                complexSquare(lMixtureMu3[currlIdx] + vThroughputMu[2]));
            complexIncrease(u + uIdx, (lMixtureAlpha[mixtureIdx] * randomPhase[uDimIdx[4]]) * (complexExponent(lMixtureC[currlIdx] + vThroughputC + sqrtMu) / sqrtMu));
            // u[uIdx] = u[uIdx] + (lMixtureAlpha[mixtureIdx] * randPhase) * (complexExponent(lMixtureC[currlIdx] + vC + sqrtMu) / sqrtMu);
        }
    }
}

int main() {
    return 0;

}
