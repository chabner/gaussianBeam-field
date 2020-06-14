#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

__constant__ double g_up; // g_up = (1 - g*g)/(2*pi)
__constant__ double g_down1; // g_down1 = 1 + g*g
__constant__ double g_down2; // g_down2 = -2*g
__constant__ double fw;
__constant__ int32_t vSize;
__constant__ int32_t lSize;

__global__ void af_ang_vl_gpu2D(double* af_ang_vl,
    const double* v1, const double* v2,
    const double* l1, const double* l2)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < (vSize * lSize))
    {
        int lIdx, vIdx;
        double cosang;

        vIdx = uIdx % vSize;
        lIdx = uIdx / vSize;

        // evalampfunc_general
        cosang = fma(v1[vIdx], l1[lIdx], v2[vIdx] * l2[lIdx]);

        if (fw == 1)
        {
            af_ang_vl[uIdx] = sqrt(g_up / fma(g_down2, cosang, g_down1));
        }
        else
        {
            af_ang_vl[uIdx] = sqrt(fw * g_up / fma(g_down2, cosang, g_down1) +
                (1 - fw) * g_up / fma(-g_down2, cosang, g_down1));
        }
    }
}

int main() {
    return 0;
}