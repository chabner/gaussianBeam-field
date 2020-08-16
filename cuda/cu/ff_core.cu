#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238
#define MAX_DIM 32
#define MAX_IN 10

// Geometry
__constant__ double box_min[3];
__constant__ double box_max[3];
__constant__ double sigt; // sigt/2
__constant__ uint32_t dims; // 2 - 2D; 3 - 3D

// scattering data
__constant__ uint32_t scatteringType; // 1 - isotropic; 2 - tabulated; 3 - HG

// - HG
__constant__ double g_up; // g_up = (1 - g*g)/(4*pi)
__constant__ double g_down1; // g_down1 = 1 + g*g
__constant__ double g_down2; // g_down2 = -2*g
__constant__ double fw;

// - tabulated
__constant__ uint32_t nAmpfunc;

// af_ang_vl data
__constant__ uint32_t af_ang_vl_maxDim; // max dim is defined from 1 to MAX_DIM of dimProd
__constant__ uint32_t af_ang_vl_dimProd[MAX_DIM];
__constant__ uint32_t af_ang_vl_inDimProd[MAX_IN * MAX_DIM]; // product of all input data

// el data
__constant__ uint32_t el_maxDim;
__constant__ uint32_t el_dimProd[MAX_DIM];
__constant__ uint32_t el_inDimProd[MAX_IN * MAX_DIM];

// ff scattering
__constant__ uint32_t ff_scattering_maxDim;
__constant__ uint32_t ff_scattering_dimProd[MAX_DIM];
__constant__ uint32_t ff_scattering_inDimProd[MAX_IN * MAX_DIM];

// refocus
__constant__ uint32_t refocus_maxDim;
__constant__ uint32_t refocus_dimProd[MAX_DIM];
__constant__ uint32_t refocus_inDimProd[MAX_IN * MAX_DIM];
__constant__ uint32_t refocus_ffElements;

// 1 for field and 2 for correlation
__constant__ uint32_t ff_corrParamNum;

__inline __device__ double2 operator +(double2 a, double2 b)
{
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    return a;
}

__inline __device__ double2 operator +(double2 a, double b)
{
    a.x = a.x + b;
    return a;
}

__inline __device__ double2 operator /(double2 a, double2 b)
{
    double2 c;
    double denominator = 1 / (fma(b.x, b.x, b.y * b.y));
    c.x = (fma(a.x, b.x, a.y * b.y)) * denominator;
    c.y = (fma(a.y, b.x, -a.x * b.y)) * denominator;
    return c;
}

__inline __device__ double2 operator /(double a, double2 b)
{
    double2 c;
    double denominator = 1 / (fma(b.x, b.x, b.y * b.y));
    c.x = (a * b.x) * denominator;
    c.y = (-a * b.y) * denominator;
    return c;
}

__inline __device__ double2 operator *(double a, double2 b)
{
    b.x = b.x * a;
    b.y = b.y * a;
    return b;
}

__inline __device__ double2 operator *(double2 a, double2 b)
{
    double2 tmpVar = a;

    a.x = fma(a.x, b.x, -a.y * b.y);
    a.y = fma(tmpVar.x, b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ double2 conjMult(double2 a, double2 b)
{
    double2 tmpVar = a;

    a.x = fma(a.x, b.x, a.y * b.y);
    a.y = fma(tmpVar.x, -b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ double2 complexConj(double2 a)
{
    a.y = -1 * a.y;
    return a;
}

__inline __device__ double2 complexSquare(double2 a)
{
    double tmpVar = a.x;
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

__device__ double2 complexExponent(double2 a)
{
    double expr, sina, cosa;
    expr = exp(a.x);
    sincos(a.y, &sina, &cosa);
    a.x = expr * cosa;
    a.y = expr * sina;
    return a;
}

__inline __device__ void complexIncrease(double2* address, double2 val)
{
    atomicAdd(&(address->x), val.x);
    atomicAdd(&(address->y), val.y);
}

__device__ double2 complexLog(double2 val)
{
    double2 logRes;
    logRes.x = log(sqrt(val.x * val.x + val.y * val.y));
    logRes.y = atan(val.y / val.x);
    return logRes;
}

// Based on MATLAB implementation
__device__ void ind2sub(unsigned int* outInd, unsigned int* siz, int sizDim, unsigned int ndx)
{
    unsigned int vi, idx;
    // siz is already cumprod of the desired size
    for (idx = (sizDim - 1); idx >= 1; idx--)
    {
        vi = ndx % siz[idx - 1];
        outInd[idx] = (ndx - vi) / siz[idx - 1];
        ndx = vi;
    }

    outInd[0] = ndx % siz[0];
}

__device__ void sub2ind(unsigned int* siz, unsigned int sizDim, unsigned int* outVec, unsigned int outDim, unsigned int* inDimProd)
{
    unsigned int idx = 0, outNum;

    for (; idx < sizDim; idx++)
    {
        for (outNum = 0; outNum < outDim; outNum++)
        {
            outVec[outNum] += (*(inDimProd + idx + MAX_DIM * outNum)) * siz[idx];
        }
    }
}

__device__ double2 evalphaseatt(double3 v, double3 dirv, const double* x, int signz, double wavelength)
{
    double bd1, bd2, bd3, dz;
    double2 e;

    // Cube dist
    if (dirv.x <= 0)
    {
        bd1 = (x[0] - box_min[0]) / abs(dirv.x);
    }
    else
    {
        bd1 = (box_max[0] - x[0]) / abs(dirv.x);
    }

    if (dirv.y <= 0)
    {
        bd2 = (x[1] - box_min[1]) / abs(dirv.y);
    }
    else
    {
        bd2 = (box_max[1] - x[1]) / abs(dirv.y);
    }

    if (dims == 3)
    {
        if (dirv.z <= 0)
        {
            bd3 = (x[2] - box_min[2]) / abs(dirv.z);
        }
        else
        {
            bd3 = (box_max[2] - x[2]) / abs(dirv.z);
        }
    }

    dz = abs(fmin(fmin(bd1, bd2), bd3));

    e.x = -sigt * dz;
    e.y = -2 * PI / wavelength * signz * (x[0] * v.x + x[1] * v.y + x[2] * v.z);

    return complexExponent(e);
}


__inline __device__ double3 operator -(double3 a)
{
    a.x = -a.x; a.y = -a.y; a.z = -a.z;
    return a;
}

__device__ double scattering(double cosang, const double* ampfunc)
{
    if (dims == 3)
    {
        if (scatteringType == 1)
        {
            return 1.0;
        }
        else if (scatteringType == 2)
        {
            cosang = acos(cosang); // scat_angle
            return ampfunc[__double2int_rn(cosang * (nAmpfunc - 1) / (PI))];
        }
        else if (scatteringType == 3)
        {
            if (fw == 1)
            {
                return sqrt(g_up / pow(fma(g_down2, cosang, g_down1), 1.5));
            }
            else
            {
                return sqrt(fw * g_up / pow(fma(g_down2, cosang, g_down1), 1.5) +
                    (1 - fw) * g_up / pow(fma(-g_down2, cosang, g_down1), 1.5));
            }
        }
    }
    else
    {
        if (scatteringType == 1)
        {
            return 1.0;
        }
        else if (scatteringType == 2)
        {
            cosang = acos(cosang); // scat_angle
            return ampfunc[__double2int_rn(cosang * (nAmpfunc - 1) / (2 * PI))];
        }
        else if (scatteringType == 3)
        {
            if (fw == 1)
            {
                return sqrt(g_up / fma(g_down2, cosang, g_down1));
            }
            else
            {
                return sqrt(fw * g_up / fma(g_down2, cosang, g_down1) +
                    (1 - fw) * g_up / fma(-g_down2, cosang, g_down1));
            }
        }
    }
}

__global__ void af_ang_vl(double* af_ang_vl, const double* ampfunc,
    const double* v1, const double* v2, const double* v3,
    const double* l1, const double* l2, const double* l3)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < af_ang_vl_dimProd[af_ang_vl_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, af_ang_vl_dimProd, af_ang_vl_maxDim, uIdx);

        // inDimProd dims:
        // 0  - v1
        // 1  - v2
        // 2  - v3
        // 3  - l1
        // 4  - l2
        // 5  - l3

        unsigned int dimsVec[6] = { 0 };
        sub2ind(uDimIdx, af_ang_vl_maxDim, dimsVec, 6, af_ang_vl_inDimProd);

        double cosang;

        // evalampfunc_general
        if (dims == 3)
        {
            cosang = fma(v1[dimsVec[0]], l1[dimsVec[3]], fma(v2[dimsVec[1]], l2[dimsVec[4]], v3[dimsVec[2]] * l3[dimsVec[5]]));
        }
        else
        {
            cosang = fma(v1[dimsVec[0]], l1[dimsVec[3]], v2[dimsVec[1]] * l2[dimsVec[4]]);
        }

        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);

        af_ang_vl[uIdx] = scattering(cosang, ampfunc);
    }
}

__global__ void ff_calcEl(double2* e_l0, double2* e_l0_ms, const double* wavelength,
    const double* ampfunc, const double* x, const double* w0, const double* w0_p,
    const double* l1, const double* l2, const double* l3,
    const double* dir_l1, const double* dir_l2, const double* dir_l3)
{
    int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < el_dimProd[el_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, el_dimProd, el_maxDim, uIdx);

        // inDimProd dims:
        // 0  - l1
        // 1  - l2
        // 2  - l3
        // 3  - dir_l1
        // 4  - dir_l2
        // 5  - dir_l3
        // 6  - wavelength

        unsigned int dimsVec[7] = { 0 };
        sub2ind(uDimIdx, el_maxDim, dimsVec, 7, el_inDimProd);

        const double* currW = w0 + dims * uDimIdx[2]; // uDimIdx[2] is the path number
        const double* currX = x + dims * uDimIdx[2];

        double3 l, dirl;
        double2 e;
        l.x = l1[dimsVec[0]]; l.y = l2[dimsVec[1]];
        dirl.x = dir_l1[dimsVec[3]]; dirl.y = dir_l2[dimsVec[4]];

        double af_l, cosang;
        if (dims == 3)
        {
            l.z = l3[dimsVec[2]];
            dirl.z = dir_l3[dimsVec[5]];

            cosang = fma(currW[0], l.x, fma(currW[1], l.y, currW[2] * l.z));
        }
        else
        {
            cosang = fma(currW[0], l.x, currW[1] * l.y);
        }

        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);

        af_l = scattering(cosang, ampfunc) / w0_p[uDimIdx[2]];
        e = evalphaseatt(l, -dirl, currX, -1, wavelength[dimsVec[6]]);

        e_l0[uIdx] = e;
        e_l0_ms[uIdx] = af_l * e;
    }
}

__global__ void ff_single(double2* us, const double2* e_l0, const double* wavelength,
    const double* af_ang_vl, const double* x, const double2* randPhase,
    const double* v1, const double* v2, const double* v3,
    const double* dir_v1, const double* dir_v2, const double* dir_v3)
{
    unsigned int bIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (bIdx < ff_scattering_dimProd[ff_scattering_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, ff_scattering_dimProd, ff_scattering_maxDim, bIdx);

        // inDimProd dims:
        // 0  - v1
        // 1  - v2
        // 2  - v3
        // 3  - dir_v1
        // 4  - dir_v2
        // 5  - dir_v3
        // 6  - af_ang_vl
        // 7  - e_l0
        // 8  - wavelength

        unsigned int dimsVec[9] = { 0 };
        sub2ind(uDimIdx, ff_scattering_maxDim, dimsVec, 9, ff_scattering_inDimProd);
        unsigned int uIdx = bIdx / ff_scattering_dimProd[0];

        double3 v, dirv;
        double2 e, field_u[2];
        double curr_af_ang_vl = af_ang_vl[dimsVec[6]];

        const double* currX = x + dims * uDimIdx[0]; // uDimIdx[0] is the path number
        double2 curr_randPhase = randPhase[uDimIdx[0]];

        v.x = v1[dimsVec[0]]; v.y = v2[dimsVec[1]];
        dirv.x = dir_v1[dimsVec[3]]; dirv.y = dir_v2[dimsVec[4]];

        if (dims == 3)
        {
            v.z = v3[dimsVec[2]];
            dirv.z = dir_v3[dimsVec[5]];
        }

        e = evalphaseatt(v, dirv, currX, 1, wavelength[dimsVec[8]]);
        field_u[0] = curr_af_ang_vl * e * curr_randPhase * e_l0[dimsVec[7]];

        if (ff_corrParamNum == 2)
        {
            v.x = v1[dimsVec[0] + 1]; v.y = v2[dimsVec[1] + 1];
            dirv.x = dir_v1[dimsVec[3] + 1]; dirv.y = dir_v2[dimsVec[4] + 1];

            if (dims == 3)
            {
                v.z = v3[dimsVec[2] + 1];
                dirv.z = dir_v3[dimsVec[5] + 1];
            }

            e = evalphaseatt(v, dirv, currX, 1, wavelength[dimsVec[8]]);

            field_u[1] = curr_af_ang_vl * e * curr_randPhase * e_l0[dimsVec[7] + 1];
            complexIncrease(us + uIdx, conjMult(field_u[0], field_u[1]));
        }
        else
        {
            complexIncrease(us + uIdx, field_u[0]);
        }

    }
}

__global__ void ff_multiple(double2* u, const double2* e_l0_ms, const double* wavelength, const bool* activatedPaths,
    const double* ampfunc, const double* x, const double* w, const double2* randPhase,
    const double* v1, const double* v2, const double* v3,
    const double* dir_v1, const double* dir_v2, const double* dir_v3)
{
    unsigned int bIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (bIdx < ff_scattering_dimProd[ff_scattering_maxDim - 1] && activatedPaths[bIdx % ff_scattering_dimProd[0]])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, ff_scattering_dimProd, ff_scattering_maxDim, bIdx);

        // inDimProd dims:
        // 0  - v1
        // 1  - v2
        // 2  - v3
        // 3  - dir_v1
        // 4  - dir_v2
        // 5  - dir_v3
        // 7  - e_l0_ms
        // 8  - wavelength

        unsigned int dimsVec[9] = { 0 };
        sub2ind(uDimIdx, ff_scattering_maxDim, dimsVec, 9, ff_scattering_inDimProd);
        unsigned int uIdx = bIdx / ff_scattering_dimProd[0];
        double3 v, dirv;
        double af_l, cosang;
        double2 e, field_u[2];

        const double* currW = w + dims * uDimIdx[0]; // uDimIdx[0] is the path number
        const double* currX = x + dims * uDimIdx[0];

        double2 curr_randPhase = randPhase[uDimIdx[0]];

        v.x = v1[dimsVec[0]]; v.y = v2[dimsVec[1]];
        dirv.x = dir_v1[dimsVec[3]]; dirv.y = dir_v2[dimsVec[4]];

        if (dims == 3)
        {
            v.z = v3[dimsVec[2]];
            dirv.z = dir_v3[dimsVec[5]];

            cosang = fma(currW[0], v.x, fma(currW[1], v.y, currW[2] * v.z));
        }
        else
        {
            cosang = fma(currW[0], v.x, currW[1] * v.y);
        }

        cosang = fmin(cosang, 0.99999999999);
        cosang = fmax(cosang, -0.99999999999);

        af_l = scattering(cosang, ampfunc);
        e = evalphaseatt(v, dirv, currX, 1, wavelength[dimsVec[8]]);
        field_u[0] = af_l * e * e_l0_ms[dimsVec[7]] * curr_randPhase;

        if (ff_corrParamNum == 2)
        {
            v.x = v1[dimsVec[0] + 1]; v.y = v2[dimsVec[1] + 1];
            dirv.x = dir_v1[dimsVec[3] + 1]; dirv.y = dir_v2[dimsVec[4] + 1];

            if (dims == 3)
            {
                v.z = v3[dimsVec[2] + 1];
                dirv.z = dir_v3[dimsVec[5] + 1];

                cosang = fma(currW[0], v.x, fma(currW[1], v.y, currW[2] * v.z));
            }
            else
            {
                cosang = fma(currW[0], v.x, currW[1] * v.y);
            }

            cosang = fmin(cosang, 0.99999999999);
            cosang = fmax(cosang, -0.99999999999);

            af_l = scattering(cosang, ampfunc);
            e = evalphaseatt(v, dirv, currX, 1, wavelength[dimsVec[8]]);
            field_u[1] = af_l * e * e_l0_ms[dimsVec[7] + 1] * curr_randPhase;

            complexIncrease(u + uIdx, conjMult(field_u[0], field_u[1]));
        }
        else
        {
            complexIncrease(u + uIdx, field_u[0]);
        }
    }
}

__global__ void refocus_correlation(double2* u_nf, const double2* u_wl, const double2* w_v)
{
    unsigned int uIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (uIdx < refocus_dimProd[refocus_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, refocus_dimProd, refocus_maxDim, uIdx);

        // inDimProd dims:
        // 0  - w_v
        // 1  - u_wl
        // 2  - w_v correlation
        // 3  - u_wl correlation

        unsigned int dimsVec[4] = { 0 };
        sub2ind(uDimIdx, refocus_maxDim, dimsVec, 2, refocus_inDimProd);
        dimsVec[2] = dimsVec[0] + refocus_ffElements;
        dimsVec[3] = dimsVec[1] + refocus_ffElements;

        double2 field_u[2];
        double2 l1, l2, v1, v2;

        field_u[0].x = 0; field_u[0].y = 0; field_u[1] = field_u[0];

        for (unsigned int idx = 0; idx < refocus_ffElements; idx++)
        {
            v1 = w_v[dimsVec[0]];
            l1 = u_wl[dimsVec[1]];
            v2 = w_v[dimsVec[2]];
            l2 = u_wl[dimsVec[3]];

            field_u[0] = v1 * l1 + field_u[0];
            field_u[1] = v2 * l2 + field_u[1];

            dimsVec[0]++;
            dimsVec[1]++;
            dimsVec[2]++;
            dimsVec[3]++;
        }

        u_nf[uIdx] = u_nf[uIdx] + conjMult(field_u[0], field_u[1]);
    }
}

__global__ void refocus_field(double2* u_nf, const double2* u_wl, const double2* w_v)
{
    unsigned int uIdx = blockIdx.x * blockDim.x + threadIdx.x;

    if (uIdx < refocus_dimProd[refocus_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, refocus_dimProd, refocus_maxDim, uIdx);

        // inDimProd dims:
        // 0  - w_v
        // 1  - u_wl

        unsigned int dimsVec[2] = { 0 };
        sub2ind(uDimIdx, refocus_maxDim, dimsVec, 2, refocus_inDimProd);

        double2 field_u;
        double2 l, v;

        field_u.x = 0; field_u.y = 0;

        for (unsigned int idx = 0; idx < refocus_ffElements; idx++)
        {
            v = w_v[dimsVec[0]];
            l = u_wl[dimsVec[1]];

            field_u = v * l + field_u;

            dimsVec[0]++;
            dimsVec[1]++;
        }

        u_nf[uIdx] = u_nf[uIdx] + field_u;
    }
}

int main() {
    return 0;

}
