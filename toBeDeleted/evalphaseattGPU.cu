#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

__constant__ double box_min[3];
__constant__ double box_max[3];
__constant__ double x[3];
__constant__ double sigt; // sigt/2
__constant__ int32_t vSize;

__global__ void evalphaseattGPU(double2* e_v0,
    const double* v1, const double* v2, const double* v3,
    const double* dir_v1, const double* dir_v2, const double* dir_v3)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < vSize)
    {
        double bd1, bd2, bd3, d;
        double v_val_1 = v1[uIdx], v_val_2 = v2[uIdx], v_val_3 = v3[uIdx];
        double v_dir_1 = dir_v1[uIdx], v_dir_2 = dir_v2[uIdx], v_dir_3 = dir_v3[uIdx];
        double att, phase, sinptr, cosptr;

        // Cube dist
        if (v_dir_1 <= 0)
        {
            bd1 = (x[0] - box_min[0]) / abs(v_dir_1);
        }
        else
        {
            bd1 = (box_max[0] - x[0]) / abs(v_dir_1);
        }

        if (v_dir_2 <= 0)
        {
            bd2 = (x[1] - box_min[1]) / abs(v_dir_2);
        }
        else
        {
            bd2 = (box_max[1] - x[1]) / abs(v_dir_2);
        }

        if (v_dir_3 <= 0)
        {
            bd3 = (x[2] - box_min[2]) / abs(v_dir_3);
        }
        else
        {
            bd3 = (box_max[2] - x[2]) / abs(v_dir_3);
        }

        d = fmin(fmin(bd1, bd2), bd3);

        // evalphaseatt
        phase = -2 * fma(v_val_1, x[0], fma(v_val_2, x[1], v_val_3 * x[2]));
        att = exp(-sigt * d);
        sincospi(phase, &sinptr, &cosptr);

        e_v0[uIdx].x = att * cosptr;
        e_v0[uIdx].y = att * sinptr;
    }
}

int main() {
    return 0;
}