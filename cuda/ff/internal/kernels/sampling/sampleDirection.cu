#include "ff_sampling.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

template<typename T, typename T3>
__global__ void sampleDirection_hg(T3 *w0, T *pw0, const importanceSampling<T,T3>* iS, ub32 dims, ub64 seed)
{
    T3 w, dir_l;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    // Choose a random direction among possible directions
    ub32 current_n = curand(&state) % iS->directions_num;
    w = *(iS->directions_list + current_n);
    dir_l = w;

    if(dims == 3)
    {
        sample_direction<T,T3,3>(&w, iS->g0_sampling, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,3,false>(iS->g0_sampling, &w, &dir_l);
    }
    if(dims == 2)
    {
        sample_direction<T,T3,2>(&w, iS->g0_sampling, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,2,false>(iS->g0_sampling, &w, &dir_l);
    }
    
    w0[threadIdx.x] = w;
}

template<typename T, typename T3>
__global__ void sampleDirection_tabulated(T3 *w0, T *pw0, const importanceSampling<T,T3>* iS, ub32 dims, ub64 seed)
{
    T3 w, dir_l;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    // Choose a random direction among possible directions
    ub32 current_n = curand(&state) % iS->directions_num;
    w = *(iS->directions_list + current_n);
    dir_l = w;

    if(dims == 3)
    {
        sample_direction<T,T3,3>(&w, iS->f0_sampling, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,3,false>(iS->f0_sampling, &w, &dir_l);
    }
    if(dims == 2)
    {
        sample_direction<T,T3,2>(&w, iS->f0_sampling, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,2,false>(iS->f0_sampling, &w, &dir_l);
    }
    
    w0[threadIdx.x] = w;
}

template<typename T, typename T3>
__global__ void sampleDirection_random(T3 *w0, T *pw0, ub32 dims, ub64 seed)
{
    T3 w;
    randomScatter scattering_structe = 0;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    if(dims == 3)
    {
        sample_direction<T,T3,3>(&w, &scattering_structe, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,3,false>(&scattering_structe, &w, &w);
    }
    if(dims == 2)
    {
        sample_direction<T,T3,2>(&w, &scattering_structe, &state);
        pw0[threadIdx.x] = scattering_contribution<T,T3,2,false>(&scattering_structe, &w, &w);
    }

    w0[threadIdx.x] = w;
}

template<typename T, typename T3>
void sampleDirection(T3 *w0, T* pw0,
        const importanceSampling<T,T3>* iS,
        ub32 is_correlation, ub32 sample_direction_flag, ub32 dims,
        const void* globalMem, ub64 curr_seed)
{
    if(sample_direction_flag == 1)
    {
        // uniform sampling
        sampleDirection_random<T,T3><<<1, THREADS_NUM>>>(w0, pw0, dims, curr_seed);
    }
    else if(sample_direction_flag == 2)
    {
        // tabulated sampling
        sampleDirection_tabulated<T,T3><<<1, THREADS_NUM>>>(w0, pw0, iS, dims, curr_seed);
    }
    else if(sample_direction_flag == 3)
    {
        // HG g0 sampling
        sampleDirection_hg<T,T3><<<1, THREADS_NUM>>>(w0, pw0, iS, dims, curr_seed);
    }

    if(!is_correlation)
    {
        // in field we take the sqrt of pw0
        global_sqrt<T><<<1, THREADS_NUM>>>(pw0);
    }
}

template void sampleDirection<double,double3>(double3 *w0, double* pw0, const importanceSampling<double,double3>* iS,
        ub32 is_correlation, ub32 sample_direction_flag, ub32 dims, const void* globalMem, ub64 curr_seed);

template void sampleDirection<float,float3>(float3 *w0, float* pw0, const importanceSampling<float,float3>* iS,
        ub32 is_correlation, ub32 sample_direction_flag, ub32 dims, const void* globalMem, ub64 curr_seed);
