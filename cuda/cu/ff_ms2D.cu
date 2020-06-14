#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

__constant__ double box_min[2];
__constant__ double box_max[2];
// __constant__ double x[2];
// __constant__ double w[2];
__constant__ double g_up; // g_up = (1 - g*g)/(2*pi)
__constant__ double g_down1; // g_down1 = 1 + g*g
__constant__ double g_down2; // g_down2 = -2*g
__constant__ double fw;
__constant__ double sigt; // sigt/2
// __constant__ double2 constCont; // sqrt(weight./px) * exp(2*pi*1i*rand)
__constant__ int32_t vSize;
__constant__ int32_t lSize;

__constant__ double fastConstCopy[6];

__global__ void ff_ms2D(double2* u, double2* e_l0_ms,
    const double* v1, const double* v2,
    const double* dir_v1, const double* dir_v2)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < (vSize * lSize))
    {
        int lIdx, vIdx;
        double bd1, bd2, d;
        double v_val_1, v_val_2;
        double v_dir_1, v_dir_2;
        double att, phase, sinptr, cosptr, cosang;
        double af_v;
        double2 l_mult_const;
        double* x, * w, * constCont;

        x = fastConstCopy;
        w = fastConstCopy + 2;
        constCont = fastConstCopy + 4;

        vIdx = uIdx % vSize;
        lIdx = uIdx / vSize;

        v_val_1 = v1[vIdx];
        v_val_2 = v2[vIdx];

        v_dir_1 = dir_v1[vIdx];
        v_dir_2 = dir_v2[vIdx];

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

        d = fmin(bd1, bd2);

        // evalphaseatt
        phase = -2 * fma(v_val_1, x[0], v_val_2 * x[1]);
        att = exp(-sigt * d);
        sincospi(phase, &sinptr, &cosptr);

        // evalampfunc_general
        cosang = fma(v_val_1, w[0], v_val_2 * w[1]);

        if (fw == 1)
        {
            af_v = sqrt(g_up / fma(g_down2, cosang, g_down1));
        }
        else
        {
            af_v = sqrt(fw * g_up / fma(g_down2, cosang, g_down1) +
                (1 - fw) * g_up / fma(-g_down2, cosang, g_down1));
        }

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