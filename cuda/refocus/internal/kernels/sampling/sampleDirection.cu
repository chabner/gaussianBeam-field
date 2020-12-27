#include "refocus_sampling.h"

template<typename T, typename T3>
__global__ void sampleDirection_random(T3 *w0, T *pw0, ub64 seed)
{
    randomScatter scattering_structe = 0;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    sample_direction<T,T3,3>(w0 + threadIdx.x, &scattering_structe, &state);
    pw0[threadIdx.x] = scattering_contribution<T,T3,3,false>(&scattering_structe, w0, w0);
}

template<typename T, typename T3>
void sampleDirection(T3 *w0, T* pw0, ub32 is_correlation, ub64 curr_seed)
{
    sampleDirection_random<T,T3><<<1, THREADS_NUM>>>(w0,pw0,curr_seed);

    if(!is_correlation)
    {
        // in field we take the sqrt of pw0
        global_sqrt<T><<<1, THREADS_NUM>>>(pw0);
    }
}

template void sampleDirection<double,double3>(double3 *w0, double* pw0, ub32 is_correlation, ub64 curr_seed);
template void sampleDirection<float,float3>(float3 *w0, float* pw0, ub32 is_correlation, ub64 curr_seed);
