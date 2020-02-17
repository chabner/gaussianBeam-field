#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <math.h>
#define PI 3.14159265358979323846

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
    c.x = (a.x * b.x + a.y * b.y) / (b.x * b.x + b.y * b.y);
    c.y = (a.y * b.x - a.x * b.y) / (b.x * b.x + b.y * b.y);
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

    a.x = a.x * b.x - a.y * b.y;
    a.y = tmpVar.x * b.y + tmpVar.y * b.x;
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
    r = sqrt(a.x * a.x + a.y * a.y);
    absrz = sqrt(r) / (sqrt((r + a.x) * (r + a.x) + a.y * a.y));
    a.x = absrz * (a.x + r);
    a.y = absrz * a.y;
    return a;
}

__forceinline __device__ double2 complexExponent(double2 a)
{
    register double expr;
    expr = exp(a.x);
    a.x = expr * cos(a.y);
    a.y = expr * sin(a.y);
    return a;
}

__forceinline __device__ double2 complexConj(double2 a)
{
    a.y = -a.y;
    return a;
}

__global__ void integrateMult( double2 * u, const int* uDim, const double2* randPhase, const double2* el, const int* elDim,
    const double2* thVMu1, const double2* thVMu2, const double2* thVMu3, const double2* thVC, const int* thVDim ,
    const double2* mixtureMu1, const double2* mixtureMu2, const double2* mixtureMu3, const double2* mixtureC, const double* mixtureAlpha, const int* mixtureDim)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx >= uDim[0] * uDim[1] * uDim[2])
        return;
    int3 uDimIdx;
    uDimIdx.z = uIdx / (uDim[0] * uDim[1]);
    uDimIdx.y = (uIdx - uDimIdx.z * (uDim[0] * uDim[1])) / uDim[0];
    uDimIdx.x = uIdx % uDim[0];
    int thVIdx = (thVDim[1] > 1) * uDimIdx.x + thVDim[1] * (thVDim[2] > 1)* uDimIdx.y + thVDim[1] * thVDim[2] * (thVDim[3] > 1)* uDimIdx.z;
    int eIdx = (elDim[1] > 1)* uDimIdx.x + elDim[1] * (elDim[2] > 1)* uDimIdx.y + elDim[1] * elDim[2] * (elDim[3] > 1)* uDimIdx.z;
    double2 sqrtMu;
    
    for (int mixtureIdx = 0; mixtureIdx < *mixtureDim; mixtureIdx++)
    {
        sqrtMu = complexSqrt(complexSquare(thVMu1[thVIdx] + complexConj(mixtureMu1[mixtureIdx])) +
                             complexSquare(thVMu2[thVIdx] + complexConj(mixtureMu2[mixtureIdx])) +
                             complexSquare(thVMu3[thVIdx] + complexConj(mixtureMu3[mixtureIdx])));
        u[uIdx] = u[uIdx] + (mixtureAlpha[mixtureIdx] * 2 * PI * *randPhase) * (el[eIdx] * complexExponent(thVC[thVIdx] + complexConj(mixtureC[mixtureIdx]) + sqrtMu) / sqrtMu);
    }

}

int main(){
    return 0;

}