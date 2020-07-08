#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "stdint.h"

#define PI 3.141592653589793238
#define MAX_DIM 32
#define MAX_IN 15

// conv specific data
__constant__ uint32_t conv_maxDim; // max dim is defined from 1 to MAX_DIM of dimProd
__constant__ uint32_t conv_dimProd[MAX_DIM];
__constant__ uint32_t conv_inDimProd[MAX_IN * MAX_DIM]; // product of all input data

// el specific data
__constant__ uint32_t el_dimProd[MAX_DIM];
__constant__ uint32_t el_inDimProd[MAX_IN * MAX_DIM];

// single scattering specific data
__constant__ uint32_t scatt_maxDim;
__constant__ uint32_t scatt_dimProd[MAX_DIM];
__constant__ uint32_t scatt_inDimProd[MAX_IN * MAX_DIM];

// multiple single scattering specific data
__constant__ uint32_t mscatt_inDimProd[MAX_IN * MAX_DIM];

// Geometry
__constant__ double box_min[3];
__constant__ double box_max[3];
__constant__ double sigt; // sigt/2

// vMF mixture data
__constant__ double mixtureAlpha[MAX_DIM];
__constant__ double mixtureMu[MAX_DIM];
__constant__ double mixtureC[MAX_DIM];
__constant__ uint32_t mixturesNum;

// 1 for field and 2 for correlation
__constant__ uint32_t corrParamNum;


__inline __device__ double2 operator +(double2 a, double2 b)
{
    a.x = a.x + b.x;
    a.y = a.y + b.y;
    return a;
}

__inline __device__ double2 operator -(double2 a, double2 b)
{
    a.x = a.x - b.x;
    a.y = a.y - b.y;
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

__inline __device__ double2 operator /(double2 a, double b)
{
    a.x = a.x / b;
    a.y = a.y / b;
    return a;
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

__inline __device__ double2 wrapTo2Pi(double2 a)
{
    a.y = fmod(a.y,2*PI);
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
    logRes.y = atan2(val.y , val.x);
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

__device__ void apToth(double2* mu, double3 v, const double* x, int signz)
{
    double bd1, bd2, bd3, dz;

    // Cube dist
    if (v.x <= 0)
    {
        bd1 = (x[0] - box_min[0]) / abs(v.x);
    }
    else
    {
        bd1 = (box_max[0] - x[0]) / abs(v.x);
    }

    if (v.y <= 0)
    {
        bd2 = (x[1] - box_min[1]) / abs(v.y);
    }
    else
    {
        bd2 = (box_max[1] - x[1]) / abs(v.y);
    }

    if (v.z <= 0)
    {
        bd3 = (x[2] - box_min[2]) / abs(v.z);
    }
    else
    {
        bd3 = (box_max[2] - x[2]) / abs(v.z);
    }


    dz = abs(fmin(fmin(bd1, bd2), bd3));

    // Throughput
    mu[0].y -= 2 * signz * PI * x[0];
    mu[1].y -= 2 * signz * PI * x[1];
    mu[2].y -= 2 * signz * PI * x[2];
    mu[3] = mu[3] + (-sigt * dz);
}


__inline __device__ double3 operator -(double3 a)
{
    a.x = -a.x; a.y = -a.y; a.z = -a.z;
    return a;
}

__device__ double2 calcSqrtMu(double2 a, double2 b, double2 c)
{
    double x_m, y_m;

    x_m = (a.x + a.y) * (a.x - a.y) + (b.x + b.y) * (b.x - b.y) + (c.x + c.y) * (c.x - c.y);
    y_m = 2.0 * (a.x * a.y + b.x * b.y + c.x * c.y);

    // r is b.x
    // 1/r is b.y
    // absrz is c.x

    b.x = sqrt(fma(x_m , x_m , y_m * y_m));
    b.y = 1 / b.x;

    c.x = 1 / sqrt(( fma(b.x + x_m , b.x + x_m , y_m * y_m )) * b.y);
    a.x = c.x * (x_m + b.x);
    a.y = c.x * y_m;

    return a;
}

__global__ void movmfConv(double2* mu1, double2* mu2, double2* mu3, double2* c,
    double2* lThMu1, double2* lThMu2, double2* lThMu3, double2* lThC, const double* x,
    const double2* lApMu1, const double2* lApMu2, const double2* lApMu3, const double2* lApC,
    const double* dirv_x, const double* dirv_y, const double* dirv_z,
    const double* dirl_x, const double* dirl_y, const double* dirl_z)
{
    unsigned int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < conv_dimProd[conv_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, conv_dimProd, conv_maxDim, uIdx);

        // inDimProd dims:
        // 0  - lApMu1
        // 1  - lApMu2
        // 2  - lApMu3
        // 3  - lApC
        // 4  - dirv_x
        // 5  - dirv_y
        // 6  - dirv_z
        // 7  - lThMu1
        // 8  - lThMu2
        // 9  - lThMu3
        // 10 - lThC
        // 11 - dirl_x
        // 12 - dirl_y
        // 13 - dirl_z

        unsigned int dimsVec[14] = { 0 };
        sub2ind(uDimIdx, conv_maxDim, dimsVec, 14, conv_inDimProd);

        double3 v_l, v_v, w;
        double2 v_mu_c[4];
        const double* currX = x + uDimIdx[2] * 3;

        v_v.x = dirv_x[dimsVec[4]]; v_v.y = dirv_y[dimsVec[5]]; v_v.z = dirv_z[dimsVec[6]];
        v_l.x = dirl_x[dimsVec[11]]; v_l.y = dirl_y[dimsVec[12]]; v_l.z = dirl_z[dimsVec[13]];

        v_mu_c[0] = lApMu1[dimsVec[0]]; v_mu_c[1] = lApMu2[dimsVec[1]]; v_mu_c[2] = lApMu3[dimsVec[2]]; v_mu_c[3] = lApC[dimsVec[3]];
        apToth(v_mu_c, -v_l, currX, -1);

        double gamma_s = mixtureMu[uDimIdx[0]];

        double2 beta_0 = complexSqrt(
            complexSquare(v_mu_c[0] + gamma_s * v_v.x) +
            complexSquare(v_mu_c[1] + gamma_s * v_v.y) +
            complexSquare(v_mu_c[2] + gamma_s * v_v.z));

        double2 gamma_s_over_beta_0 = gamma_s / beta_0;

        mu1[uIdx] = gamma_s_over_beta_0 * v_mu_c[0];
        mu2[uIdx] = gamma_s_over_beta_0 * v_mu_c[1];
        mu3[uIdx] = gamma_s_over_beta_0 * v_mu_c[2];

        double real_abs_mu;

        if (abs(gamma_s) > 0.0000001)
        {
            real_abs_mu = rnorm3d(mu1[uIdx].x, mu2[uIdx].x, mu3[uIdx].x);
            w.x = mu1[uIdx].x * real_abs_mu; w.y = mu2[uIdx].x * real_abs_mu; w.z = mu3[uIdx].x * real_abs_mu;
        }
        else
        {
            w.x = 0.0; w.y = 0.0; w.z = 1.0;
        }

        double2 log_estimated_conv_max = w.x * mu1[uIdx] + w.y * mu2[uIdx] + w.z * mu3[uIdx];
        double2 C = calcSqrtMu(
            v_mu_c[0] + gamma_s * w.x,
            v_mu_c[1] + gamma_s * w.y,
            v_mu_c[2] + gamma_s * w.z
            );

        c[uIdx] = v_mu_c[3] + C + log(2 * PI) - complexLog(C) + mixtureC[uDimIdx[0]] - log_estimated_conv_max;
        lThMu1[dimsVec[7]] = v_mu_c[0];
        lThMu2[dimsVec[8]] = v_mu_c[1];
        lThMu3[dimsVec[9]] = v_mu_c[2];
        lThC[dimsVec[10]] = v_mu_c[3];
    }
}

__global__ void integrateEl(double2* el, const double* w, const double* w0p,
    const double2* lThMu1, const double2* lThMu2, const double2* lThMu3, const double2* lThC)
{
    unsigned int uIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (uIdx < el_dimProd[conv_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, el_dimProd, conv_maxDim, uIdx);

        // inDimProd dims:
        // 0  - lThMu1
        // 1  - lThMu2
        // 2  - lThMu3
        // 3  - lThC

        unsigned int dimsVec[4] = { 0 }, mixtureIdx;
        sub2ind(uDimIdx, conv_maxDim, dimsVec, 4, el_inDimProd);
        double2 sqrtMu;
        const double* currW = w + uDimIdx[2] * 3;
        double logw0 = log(w0p[uDimIdx[2]]);

        for (mixtureIdx = 0; mixtureIdx < mixturesNum; mixtureIdx++)
        {
            sqrtMu = calcSqrtMu(
                lThMu1[dimsVec[0]] + mixtureMu[mixtureIdx] * currW[0],
                lThMu2[dimsVec[1]] + mixtureMu[mixtureIdx] * currW[1],
                lThMu3[dimsVec[2]] + mixtureMu[mixtureIdx] * currW[2]);
            el[uIdx] = el[uIdx] + 2 * PI * mixtureAlpha[mixtureIdx] * (complexExponent(lThC[dimsVec[3]] + (mixtureC[mixtureIdx] - logw0) + sqrtMu - complexLog(sqrtMu)));
        }
        el[uIdx] = complexLog(el[uIdx]);
    }
}

__global__ void singleScattering(double2* u, const double* x, const double2* randomPhase,
    const double2* lThMu1, const double2* lThMu2, const double2* lThMu3, const double2* lThC,
    const double2* vApMu1, const double2* vApMu2, const double2* vApMu3, const double2* vApC,
    const double* dirv_x, const double* dirv_y, const double* dirv_z)
{
    unsigned int bIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (bIdx < scatt_dimProd[scatt_maxDim - 1])
    {
        unsigned int uDimIdx[MAX_DIM], mixtureIdx, lIdx, corrNum, corrIdx;
        ind2sub(uDimIdx, scatt_dimProd, scatt_maxDim, bIdx);

        // inDimProd dims:
        // 0 - lIdx
        // 1 - vIdxMu1
        // 2 - vIdxMu2
        // 3 - vIdxMu3
        // 4 - vIdxC
        // 5 - vIdxFocalDirection_x
        // 6 - vIdxFocalDirection_y
        // 7 - vIdxFocalDirection_z

        unsigned int dimsVec[8] = { 0 };
        sub2ind(uDimIdx, scatt_maxDim, dimsVec, 8, scatt_inDimProd);
        unsigned int uIdx = bIdx / scatt_dimProd[0];


        double3 v;
        double2 randPhase = randomPhase[uDimIdx[0]];

        double2 v_mu_c[8];
        double2 sqrtMu, alpha;
        double2 field_u[2];
        const double* currX = x + uDimIdx[0] * 3;

        field_u[0].x = 0; field_u[0].y = 0; field_u[1] = field_u[0];
        v.x = dirv_x[dimsVec[5]]; v.y = dirv_y[dimsVec[6]]; v.z = dirv_z[dimsVec[7]];
        v_mu_c[0] = vApMu1[dimsVec[1]]; v_mu_c[1] = vApMu2[dimsVec[2]]; v_mu_c[2] = vApMu3[dimsVec[3]]; v_mu_c[3] = vApC[dimsVec[4]];
        apToth(v_mu_c, v, currX, 1);

        if (corrParamNum == 2)
        {
            v.x = dirv_x[dimsVec[5] + 1]; v.y = dirv_y[dimsVec[6] + 1]; v.z = dirv_z[dimsVec[7] + 1];
            v_mu_c[4] = vApMu1[dimsVec[1] + 1]; v_mu_c[5] = vApMu2[dimsVec[2] + 1]; v_mu_c[6] = vApMu3[dimsVec[3] + 1]; v_mu_c[7] = vApC[dimsVec[4] + 1];
            apToth(v_mu_c + 4, v, currX, 1);
        }

        lIdx = dimsVec[0];

        for (mixtureIdx = 0; mixtureIdx < mixturesNum; mixtureIdx++)
        {
            alpha = mixtureAlpha[mixtureIdx] * randPhase;
            for (corrNum = 0; corrNum < corrParamNum; corrNum++)
            {
                corrIdx = lIdx + mixturesNum * corrNum;
                sqrtMu = calcSqrtMu(
                    lThMu1[corrIdx] + v_mu_c[0 + 4 * corrNum],
                    lThMu2[corrIdx] + v_mu_c[1 + 4 * corrNum],
                    lThMu3[corrIdx] + v_mu_c[2 + 4 * corrNum]);
                field_u[corrNum] = field_u[corrNum] + alpha * (complexExponent((lThC[corrIdx]) + (v_mu_c[3 + 4 * corrNum]) + (sqrtMu) - complexLog(sqrtMu)));
                //field_u[corrNum] = field_u[corrNum] + alpha * (complexExponent((lThC[corrIdx]) + (v_mu_c[3 + 4 * corrNum]) + (sqrtMu)));
            }
            lIdx++;
        }

        if (corrParamNum == 2)
            complexIncrease(u + uIdx, conjMult(field_u[0], field_u[1]));
        else
            complexIncrease(u + uIdx, field_u[0]);

    }
}

__global__ void multipleScattering(double2* u, const double2* el, const double* x, const double* w, const double2* randomPhase, const bool* activatedPaths,
    const double2* vApertureMu1, const double2* vApertureMu2, const double2* vApertureMu3, const double2* vApertureC,
    const double* dirv_x, const double* dirv_y, const double* dirv_z)
{
    unsigned int bIdx = blockIdx.x * blockDim.x + threadIdx.x;
    if (bIdx < scatt_dimProd[scatt_maxDim - 1] && activatedPaths[bIdx % scatt_dimProd[0]])
    {


        unsigned int uDimIdx[MAX_DIM];
        ind2sub(uDimIdx, scatt_dimProd, scatt_maxDim, bIdx);

        // inDimProd dims:
        // 0 - el
        // 1 - vApertureMu1
        // 2 - vApertureMu2
        // 3 - vApertureMu3
        // 4 - vApertureC
        // 5 - dirv_x
        // 6 - dirv_y
        // 7 - dirv_z

        unsigned int dimsVec[8] = { 0 };
        sub2ind(uDimIdx, scatt_maxDim, dimsVec, 8, mscatt_inDimProd);
        unsigned int uIdx = bIdx / scatt_dimProd[0], mixtureIdx, corrNum;

        double3 v;
        const double2* elVal = el + dimsVec[0];

        double2 randPhase = randomPhase[uDimIdx[0]], v_mu_c[8], sqrtMu, alpha, field_u[2];
        const double* currX = x + uDimIdx[0] * 3;
        const double* currW = w + uDimIdx[0] * 3;
        
        field_u[0].x = 0; field_u[0].y = 0; field_u[1] = field_u[0];

        v.x = dirv_x[dimsVec[5]]; v.y = dirv_y[dimsVec[6]]; v.z = dirv_z[dimsVec[7]];
        v_mu_c[0] = vApertureMu1[dimsVec[1]]; v_mu_c[1] = vApertureMu2[dimsVec[2]]; v_mu_c[2] = vApertureMu3[dimsVec[3]]; v_mu_c[3] = vApertureC[dimsVec[4]];
        apToth(v_mu_c, v, currX, 1);
        
        if (corrParamNum == 2)
        {
            v.x = dirv_x[dimsVec[5] + 1]; v.y = dirv_y[dimsVec[6] + 1]; v.z = dirv_z[dimsVec[7] + 1];
            v_mu_c[4] = vApertureMu1[dimsVec[1] + 1]; v_mu_c[5] = vApertureMu2[dimsVec[2] + 1]; v_mu_c[6] = vApertureMu3[dimsVec[3] + 1]; v_mu_c[7] = vApertureC[dimsVec[4] + 1];
            apToth(v_mu_c + 4, v, currX, 1);
        }

        
        for (mixtureIdx = 0; mixtureIdx < mixturesNum; mixtureIdx++)
        {
            alpha = mixtureAlpha[mixtureIdx] * randPhase;
            for (corrNum = 0; corrNum < corrParamNum; corrNum++)
            {
                sqrtMu = calcSqrtMu(
                    v_mu_c[0 + 4 * corrNum] + mixtureMu[mixtureIdx] * currW[0],
                    v_mu_c[1 + 4 * corrNum] + mixtureMu[mixtureIdx] * currW[1],
                    v_mu_c[2 + 4 * corrNum] + mixtureMu[mixtureIdx] * currW[2]);
                field_u[corrNum] = field_u[corrNum] + alpha * (complexExponent(v_mu_c[3 + 4 * corrNum] + mixtureC[mixtureIdx] + sqrtMu - complexLog(sqrtMu) + elVal[corrNum]));
                // field_u[corrNum] = field_u[corrNum] + alpha * elVal[corrNum] * (complexExponent(v_mu_c[3 + 4 * corrNum] + mixtureC[mixtureIdx] + sqrtMu)) / sqrtMu;
            }
        }

        if (corrParamNum == 2)
            complexIncrease(u + uIdx, conjMult(field_u[0], field_u[1]));
        else
            complexIncrease(u + uIdx, field_u[0]);
    }
}

int main() {
    return 0;

}
