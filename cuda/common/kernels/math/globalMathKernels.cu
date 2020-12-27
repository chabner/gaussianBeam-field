// #pragma once

#include "globalMathKernels.h"
#include "cuRandoms.cu"
#include "math.h"
#include "cuMath.cu"

template<typename T>
__global__ void global_sqrt(T* x)
{
    x[threadIdx.x] = sqrt(x[threadIdx.x]);
}

template<typename T>
__global__ void global_mult(T* x, const T* y)
{
    x[threadIdx.x] *= y[threadIdx.x];
}

template<typename T>
__global__ void global_add(T* x, const T* y)
{
    x[threadIdx.x] += y[threadIdx.x];
}

template<typename T>
__global__ void global_min(T* x, T min_val)
{
    if(x[threadIdx.x] < min_val)
    {
        x[threadIdx.x] = min_val;
    }
}

template<typename T2>
__global__ void elementwise_conj_mult(T2* x, const T2* y, ub32 n)
{
    ub32 elem_num = blockIdx.x * blockDim.x + threadIdx.x;

    if(elem_num < n)
    {
        x[elem_num] = conjMult(x[elem_num], y[elem_num]);
    }
}

template<typename T2>
__global__ void elementwise_mult(T2* x, const T2* y, ub32 n)
{
    ub32 elem_num = blockIdx.x * blockDim.x + threadIdx.x;

    if(elem_num < n)
    {
        x[elem_num] = x[elem_num] * y[elem_num];
    }
}

template<typename T, typename T3, ub32 dims>
__device__ void rotate_by_theta(T3 *w, T costheta, curandState_t *state)
{
    T sintheta = sqrt(1 - costheta * costheta), phi;
    randomUniform(state, &phi);

    if(dims == 2)
    {
        sintheta = (phi > 0.5 ? 1.0 : -1.0) * sintheta;
        T3 ow = *w;
        w->x =  costheta * ow.x + sintheta * ow.y;
        w->y = -sintheta * ow.x + costheta * ow.y;
    }
            
    if(dims == 3)
    {
        T sinphi, cosphi;

        phi *= (2.0); // * PI        
        sincospi(phi, &sinphi, &cosphi );

        if(abs(w->z) > 0.999999999)
        {
            w->x = sintheta * cosphi;
            w->y = sintheta * sinphi;
            w->z = costheta * (signbit(w->z) ? -1.0 : 1.0);
        }
        else
        {
            T sintheta_w = sqrt(1.0 - w->z * w->z);
            T3 ow = *w;
            w->x =  sintheta * (ow.x * ow.z * cosphi - ow.y * sinphi) / sintheta_w + ow.x * costheta;
            w->y =  sintheta * (ow.y * ow.z * cosphi + ow.x * sinphi) / sintheta_w + ow.y * costheta;
            w->z = -sintheta * cosphi * sintheta_w + ow.z * costheta;
        }
    }
}

template<typename T,typename T3, ub32 dims>
__global__ void normalize_dir(T3 *w)
{
    if(dims == 3)
    {
        double w_x = (double) w[threadIdx.x].x;
        double w_y = (double) w[threadIdx.x].y;
        double w_z = (double) w[threadIdx.x].z;
        double norm_factor = 1/sqrt(w_x * w_x + w_y * w_y + w_z * w_z);
        w[threadIdx.x].x = (T) (w_x * norm_factor);
        w[threadIdx.x].y = (T) (w_y * norm_factor);
        w[threadIdx.x].z = (T) (w_z * norm_factor);
    }
    
    if(dims == 2)
    {
        double w_x = (double) w[threadIdx.x].x;
        double w_y = (double) w[threadIdx.x].y;
        double norm_factor = 1/sqrt(w_x * w_x + w_y * w_y);
        w[threadIdx.x].x = (T) (w_x * norm_factor);
        w[threadIdx.x].y = (T) (w_y * norm_factor);
    }
}

template<typename T>
__device__ ub32 binary_search(T *search_array, T search_target, ub32 elements_num)
{
    ub32 m, L, R;
    L = 1; R = elements_num;

    while(L <= R)
    {
        m = (L+R)/2;

        if(search_array[m - 1] < search_target)
        {
            L = m + 1;
        }
        else
        {
            R = m - 1;
        }
    }

    if(search_array[m - 1] >= search_target)
    {   
        m--;
    }
    return m;
}

template __global__ void global_sqrt<double>(double* x);
template __global__ void global_sqrt<float>(float* x);

template __global__ void global_mult<double>(double* x, const double* y);
template __global__ void global_mult<float>(float* x, const float* y);

template __global__ void global_add<double>(double* x, const double* y);
template __global__ void global_add<float>(float* x, const float* y);

template __global__ void global_min<double>(double* x, double min_val);
template __global__ void global_min<float>(float* x, float min_val);

template __global__ void elementwise_conj_mult<double2>(double2* x, const double2* y, ub32 n);
template __global__ void elementwise_conj_mult<float2>(float2* x, const float2* y, ub32 n);

template __global__ void elementwise_mult<double2>(double2* x, const double2* y, ub32 n);
template __global__ void elementwise_mult<float2>(float2* x, const float2* y, ub32 n);

template __device__ void rotate_by_theta<double,double3,2>(double3 *w, double costheta, curandState_t *state);
template __device__ void rotate_by_theta<float,float3,2>(float3 *w, float costheta, curandState_t *state);
template __device__ void rotate_by_theta<double,double3,3>(double3 *w, double costheta, curandState_t *state);
template __device__ void rotate_by_theta<float,float3,3>(float3 *w, float costheta, curandState_t *state);

template __global__ void normalize_dir<double, double3,2>(double3 *w);
template __global__ void normalize_dir<float, float3,2>(float3 *w);
template __global__ void normalize_dir<double, double3,3>(double3 *w);
template __global__ void normalize_dir<float, float3,3>(float3 *w);

template __device__ ub32 binary_search<double>(double *search_array, double search_target, ub32 elements_num);
template __device__ ub32 binary_search<float>(float *search_array, float search_target, ub32 elements_num);
