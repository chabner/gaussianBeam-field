#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

__constant__ double box_min[2];
__constant__ double box_max[2];
// __constant__ double x[3];
__constant__ double sigt; // sigt/2
// __constant__ double2 constCont; // sqrt(weight./px) * exp(2*pi*1i*rand)
__constant__ int32_t vSize;
__constant__ int32_t lSize;
__constant__ int32_t wSize;

__constant__ double fastConstCopy[4];

__global__ void ff_single2D(double2* us, double2* wl, double2* e_l0, double2* af_ang_vl,
    const double* v1, const double* v2, const double* dir_v1, const double* dir_v2)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < (vSize * wSize))
    {
        int lIdx, vIdx, wIdx;
        double bd1, bd2, d;
        double v_val_1, v_val_2;
        double v_dir_1, v_dir_2;
        double att, phase, sinptr, cosptr;
        double* x, * constCont;
        double2 vlMult, vwMult, tpath, attRotation, currAf_ang_vl;
        double wx, wy;

        x = fastConstCopy;
        constCont = fastConstCopy + 2;

        vIdx = uIdx % vSize;
        wIdx = uIdx / vSize;

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

        // result
        for (lIdx = 0; lIdx < lSize; lIdx++)
        {
            // attRotation - const
            currAf_ang_vl = af_ang_vl[vIdx + lIdx * vSize];
            attRotation.x = att * fma(currAf_ang_vl.x, constCont[0], -currAf_ang_vl.y * constCont[1]);
            attRotation.y = att * fma(currAf_ang_vl.x, constCont[1], currAf_ang_vl.y * constCont[0]);

            // af_ang_vl .* e_v0 .* e_l0
            vlMult.x = fma(cosptr, e_l0[lIdx].x, -sinptr * e_l0[lIdx].y);
            vlMult.y = fma(cosptr, e_l0[lIdx].y, sinptr * e_l0[lIdx].x);

            // permute(Wl,[3,2,1]) .* (af_ang_vl .* e_v0 .* e_l0)
            wx = (wl + wIdx + lIdx * wSize)->x;
            wy = (wl + wIdx + lIdx * wSize)->y;
            vwMult.x = fma(vlMult.x, wx, -vlMult.y * wy);
            vwMult.y = fma(vlMult.x, wy, vlMult.y * wx);

            // multiple with const
            tpath.x = fma(vwMult.x, attRotation.x, -vwMult.y * attRotation.y);
            tpath.y = fma(vwMult.x, attRotation.y, vwMult.y * attRotation.x);

            // add to us
            us[uIdx].x += tpath.x;
            us[uIdx].y += tpath.y;
        }
    }
}

int main() {
    return 0;
}