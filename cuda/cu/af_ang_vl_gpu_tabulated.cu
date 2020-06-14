#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238

__constant__ int32_t vSize;
__constant__ int32_t lSize;
__constant__ int32_t nAmpfunc;

__global__ void af_ang_vl_gpu_tabulated(double* af_ang_vl, double* ampfunc,
    const double* v1, const double* v2, const double* v3,
    const double* l1, const double* l2, const double* l3)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < (vSize * lSize))
    {
        int lIdx, vIdx;
        double cosang, scat_angle;

        vIdx = uIdx % vSize;
        lIdx = uIdx / vSize;

        // evalampfunc_general
        cosang = fma(v1[vIdx], l1[lIdx], fma(v2[vIdx], l2[lIdx], v3[vIdx] * l3[lIdx]));
        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);
        scat_angle = acos(cosang);

        af_ang_vl[uIdx] = ampfunc[__double2int_rn(scat_angle * (nAmpfunc - 1) / (PI))];
    }
}

int main() {
    return 0;
}