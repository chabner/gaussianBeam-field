#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

__constant__ int32_t uDimProd[3];
__constant__ int32_t lDim[4];
__constant__ int32_t vDim[4];
__constant__ double2 randPhase;
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
    double denominator = 1/(fma(b.x,b.x,b.y * b.y));
    c.x = (fma(a.x,b.x,a.y * b.y)) * denominator;
    c.y = (fma(a.y,b.x,-a.x * b.y)) * denominator;
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

    a.x = fma(a.x,b.x,- a.y * b.y);
    a.y = fma(tmpVar.x,b.y,tmpVar.y * b.x);
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
    r = hypot(a.x,a.y);
    absrz = sqrt(r) * rhypot(r + a.x,a.y);
    // absrz = sqrt(r) / (fma(sqrt((r + a.x)),(r + a.x),a.y * a.y));
    a.x = absrz * (a.x + r);
    a.y = absrz * a.y;
    return a;
}

__forceinline __device__ double2 complexExponent(double2 a)
{
    register double expr,sina,cosa;
    expr = exp(a.x);
    sincos(a.y, &sina, &cosa);
    a.x = expr * cosa;
    a.y = expr * sina;
    return a;
}

__forceinline __device__ void complexIncrease(double2* address, double2 val)
{
    //atomicAdd(&(*address).x, val.x);
    //atomicAdd(&(*address).y, val.y);
    (*address).x += val.x;
    (*address).y += val.y;
}

__global__ void integrateMult(double2* u,
    const double2* lMixtureMu1, const double2* lMixtureMu2, const double2* lMixtureMu3, const double2* lMixtureC,
    const double2* vMixtureMu1, const double2* vMixtureMu2, const double2* vMixtureMu3, const double2* vMixtureC)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < uDimProd[2])
    {
        int3 uDimIdx;
        uDimIdx.z = uIdx / (uDimProd[1]);
        uDimIdx.y = (uIdx - uDimIdx.z * (uDimProd[1])) / uDimProd[0];
        uDimIdx.x = uIdx % uDimProd[0];
        int lIdx = lDim[0] * (lDim[1] > 1)* uDimIdx.x + lDim[0] * lDim[1] * (lDim[2] > 1)* uDimIdx.y + lDim[0] * lDim[1] * lDim[2] * (lDim[3] > 1)* uDimIdx.z;
        int vIdx = (vDim[1] > 1)* uDimIdx.x + vDim[1] * (vDim[2] > 1)* uDimIdx.y + vDim[1] * vDim[2] * (vDim[3] > 1)* uDimIdx.z;
        double2 sqrtMu;
        int currlIdx;
        double2 vMu1 = vMixtureMu1[vIdx];
        double2 vMu2 = vMixtureMu2[vIdx];
        double2 vMu3 = vMixtureMu3[vIdx];
        double2 vC = vMixtureC[vIdx];

        for (int mixtureIdx = 0; mixtureIdx < lDim[0]; mixtureIdx++)
        {
            currlIdx = lIdx + mixtureIdx;
            sqrtMu = complexSqrt(complexSquare(lMixtureMu1[currlIdx] + vMu1) +
                complexSquare(lMixtureMu2[currlIdx] + vMu2) +
                complexSquare(lMixtureMu3[currlIdx] + vMu3));
            complexIncrease(u + uIdx, (lMixtureAlpha[mixtureIdx] * randPhase) * (complexExponent(lMixtureC[currlIdx] + vC + sqrtMu) / sqrtMu));
            // u[uIdx] = u[uIdx] + (lMixtureAlpha[mixtureIdx] * randPhase) * (complexExponent(lMixtureC[currlIdx] + vC + sqrtMu) / sqrtMu);
        }
    }
}

int main() {
    return 0;

}