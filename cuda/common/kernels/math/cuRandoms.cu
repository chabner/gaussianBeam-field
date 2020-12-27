#pragma once

#include "cuda_runtime.h"
#include "curand.h"
#include "curand_kernel.h"
#include "math.h"

__inline__ __device__ void randomUniform(curandState_t* state, double* randNum)
{
    *randNum = curand_uniform_double(state);
}

__inline__ __device__ void randomUniform(curandState_t* state, float* randNum)
{
    *randNum = curand_uniform(state);
}

__inline__ __device__ void randomNormal2(curandState_t* state, double2* randNum)
{
    *randNum = curand_normal2_double(state);
}

__inline__ __device__ void randomNormal2(curandState_t* state, float2* randNum)
{
    *randNum = curand_normal2(state);
}

__inline__ __device__ void randomNormal2(curandState_t* state, double3* randNum)
{
    double2 tmp = curand_normal2_double(state);
    randNum->x = tmp.x;
    randNum->y = tmp.y;
}

__inline__ __device__ void randomNormal2(curandState_t* state, float3* randNum)
{
    float2 tmp = curand_normal2(state);
    randNum->x = tmp.x;
    randNum->y = tmp.y;
}

__inline__ __device__ void randomNormal(curandState_t* state, double* randNum)
{
    *randNum = curand_normal_double(state);
}

__inline__ __device__ void randomNormal(curandState_t* state, float* randNum)
{
    *randNum = curand_normal(state);
}

template<typename T, typename T3>
__device__ T3 randomDirection3(curandState_t* state)
{
    T divFactor = 0.0;
    T r_divFactor;
    T3 w;

    while(divFactor < 1e-5)
    {
        randomNormal(state,&w.x);
        randomNormal(state,&w.y);
        randomNormal(state,&w.z);
        divFactor = norm3d(w.x,w.y,w.z);
    }

    r_divFactor = 1 / divFactor;
    w.x *= r_divFactor;
    w.y *= r_divFactor;
    w.z *= r_divFactor;

    return w;
}

template<typename T, typename T2>
__device__ T2 randomDirection2(curandState_t* state)
{
    T divFactor = 0.0;
    T r_divFactor;
    T2 w;

    while(divFactor < 1e-5)
    {
        randomNormal2(state,&w);
        divFactor = hypot(w.x,w.y);
    }

    r_divFactor = 1 / divFactor;
    w.x *= r_divFactor;
    w.y *= r_divFactor;

    return w;
}

template<typename T, typename T2, typename T3>
__device__ T3 random_vMF_direction(T3 center, T kappa, curandState_t* state)
{
    T b = -kappa + hypot(kappa,(T)1.0);
    T dir_x0 = (1.0-b)/(1.0+b);
    T s_c = kappa * dir_x0 + 2.0 * log(1.0 - dir_x0 * dir_x0);

    T z, u, w_z, t;

    do
    {
        randomUniform(state,&z); // beta(1,1) is uniform
        randomUniform(state,&u);
        w_z = (1.0 - (1.0+b)*z)/(1.0 - (1.0-b)*z);
        t = kappa * w_z + 2.0 * log(1.0 - dir_x0 * w_z) - s_c;
    }
    while(t < log(u));

    T2 v = randomDirection2<T,T2>(state);

    T3 sample;
    T sqrt_w = sqrt(abs(1.0 - w_z * w_z));
    sample.x = sqrt_w * v.x;
    sample.y = sqrt_w * v.y;
    sample.z = w_z;

    // house the center
    T s = center.x * center.x + center.y * center.y;

    T3 mu_housed;
    mu_housed.x = center.x;
    mu_housed.y = center.y;
    mu_housed.z = 1.0;
    
    if(abs(s) < 1e-8)
    {
        b = 0.0;
    }
    else
    {
        mu_housed.z = (center.z <= 0.0 ? center.z - 1.0 : -s/(center.z + 1.0));
        b = 2.0 * (mu_housed.z * mu_housed.z) /(s + mu_housed.z * mu_housed.z);
        mu_housed.x = mu_housed.x / mu_housed.z;
        mu_housed.y = mu_housed.y / mu_housed.z;
        mu_housed.z = 1.0;
    }
    
    T3 sample_housed;
    sample_housed.x = (1.0 - b * mu_housed.x * mu_housed.x) * sample.x +
                      (    - b * mu_housed.x * mu_housed.y) * sample.y +
                      (    - b * mu_housed.x * mu_housed.z) * sample.z;

    sample_housed.y = (    - b * mu_housed.y * mu_housed.x) * sample.x +
                      (1.0 - b * mu_housed.y * mu_housed.y) * sample.y +
                      (    - b * mu_housed.y * mu_housed.z) * sample.z;

    sample_housed.z = (    - b * mu_housed.z * mu_housed.x) * sample.x +
                      (    - b * mu_housed.z * mu_housed.y) * sample.y +
                      (1.0 - b * mu_housed.z * mu_housed.z) * sample.z;

    return sample_housed;
}
