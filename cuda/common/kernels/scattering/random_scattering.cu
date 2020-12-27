#include "random_scattering.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const randomScatter* scattering_struct, const T3* D, const T3* W)
{
    T ret_val;

    if(dims == 3)
    {
        ret_val = 1 / (4 * M_PI);
    }

    if(dims == 2)
    {
        ret_val = 1 / (2 * M_PI);
    }

    if(is_sqrt)
    {
        ret_val = sqrt(ret_val);
    }

    return ret_val;
}

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const randomScatter* sample_function, ub64 seed, ub32 dims_num)
{
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    if(dims_num == 3)
    {
        sample_direction<T,T3,3>(w + threadIdx.x, sample_function, &state);
    }
    if(dims_num == 2)
    {
        sample_direction<T,T3,2>(w + threadIdx.x, sample_function, &state);
    }
}

template<typename T, typename T3, ub32 dims>
__device__ void sample_direction(T3 *w, const randomScatter* sample_function, curandState_t *state)
{
    if(dims == 3)
    {
        *w = randomDirection3<T,T3>(state);
    }

    if(dims == 2)
    {
        *w = randomDirection2<T,T3>(state);
    }
}


template __global__ void sample_direction<double,double3>(double3 *w, const randomScatter* sample_function, ub64 seed, ub32 dims_num);
template __global__ void sample_direction<float,float3>(float3 *w, const randomScatter* sample_function, ub64 seed, ub32 dims_num);

template __device__ void sample_direction<double,double3,2>(double3 *w, const randomScatter* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,2>(float3 *w, const randomScatter* sample_function, curandState_t *state);
template __device__ void sample_direction<double,double3,3>(double3 *w, const randomScatter* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,3>(float3 *w, const randomScatter* sample_function, curandState_t *state);

template __device__ double scattering_contribution<double,double3,2,true>(const randomScatter* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,true>(const randomScatter* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,true>(const randomScatter* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,true>(const randomScatter* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,2,false>(const randomScatter* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,false>(const randomScatter* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,false>(const randomScatter* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,false>(const randomScatter* scattering_struct, const float3* D, const float3* W);
