#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238

__constant__ double box_min[2];
__constant__ double box_max[2];
// __constant__ double x[2];
// __constant__ double w[2];
__constant__ double sigt; // sigt/2
// __constant__ double2 constCont; // sqrt(weight./px) * exp(2*pi*1i*rand)
__constant__ int32_t vSize;
__constant__ int32_t lSize;
__constant__ int32_t nAmpfunc;

__constant__ double fastConstCopy[6];

__global__ void ff_ms_tabulated2D(double2* u, double2* e_l0_ms, double2* ampfunc,
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
        double scat_angle;
        double2 af_v, complex_att;
        double2 l_mult;
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
        scat_angle = acos(cosang);
        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);
        af_v = ampfunc[__double2int_rn(scat_angle * (nAmpfunc - 1) / (2 * PI))];

        // multiple with the constant
        complex_att.x = att * fma(af_v.x, constCont[0], -af_v.y * constCont[1]);
        complex_att.y = att * fma(af_v.x, constCont[1], af_v.y * constCont[0]);

        // result
        l_mult.x = fma(e_l0_ms[lIdx].x, complex_att.x, -e_l0_ms[lIdx].y * complex_att.y);
        l_mult.y = fma(e_l0_ms[lIdx].x, complex_att.y, e_l0_ms[lIdx].y * complex_att.x);
        u[uIdx].x += fma(cosptr, l_mult.x, -sinptr * l_mult.y);
        u[uIdx].y += fma(cosptr, l_mult.y, sinptr * l_mult.x);
    }
}

int main() {
    return 0;
}