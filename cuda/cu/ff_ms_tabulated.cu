#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238

__constant__ double box_min[3];
__constant__ double box_max[3];
// __constant__ double x[3];
// __constant__ double w[3];
__constant__ double sigt; // sigt/2
// __constant__ double2 constCont; // sqrt(weight./px) * exp(2*pi*1i*rand)
__constant__ int32_t vSize;
__constant__ int32_t lSize;
__constant__ int32_t nAmpfunc;

__constant__ double fastConstCopy[8];

__global__ void ff_ms_tabulated(double2* u, double2* e_l0_ms, double* ampfunc,
    const double* v1, const double* v2, const double* v3,
    const double* dir_v1, const double* dir_v2, const double* dir_v3)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < (vSize * lSize))
    {
        int lIdx, vIdx;
        double bd1, bd2, bd3, d;
        double v_val_1, v_val_2, v_val_3;
        double v_dir_1, v_dir_2, v_dir_3;
        double att, phase, sinptr, cosptr, cosang;
        double af_v, scat_angle;
        double2 l_mult_const;
        double* x, * w, * constCont;

        x = fastConstCopy;
        w = fastConstCopy + 3;
        constCont = fastConstCopy + 6;

        vIdx = uIdx % vSize;
        lIdx = uIdx / vSize;

        v_val_1 = v1[vIdx];
        v_val_2 = v2[vIdx];
        v_val_3 = v3[vIdx];

        v_dir_1 = dir_v1[vIdx];
        v_dir_2 = dir_v2[vIdx];
        v_dir_3 = dir_v3[vIdx];

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

        // evalampfunc_general
        cosang = fma(v_val_1, w[0], fma(v_val_2, w[1], v_val_3 * w[2]));
        scat_angle = acos(cosang);
        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);
        af_v = ampfunc[__double2int_rn(scat_angle * (nAmpfunc - 1) / (PI))];

        att *= af_v;

        // result
        l_mult_const.x = fma(e_l0_ms[lIdx].x, constCont[0], -e_l0_ms[lIdx].y * constCont[1]);
        l_mult_const.y = fma(e_l0_ms[lIdx].x, constCont[1], e_l0_ms[lIdx].y * constCont[0]);
        u[uIdx].x += att * fma(cosptr, l_mult_const.x, -sinptr * l_mult_const.y);
        u[uIdx].y += att * fma(cosptr, l_mult_const.y, sinptr * l_mult_const.x);
    }
}

int main() {
    return 0;
}