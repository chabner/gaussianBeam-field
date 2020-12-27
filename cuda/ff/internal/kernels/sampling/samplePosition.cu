#include "ff_sampling.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ constStructre<double> constStr_smp_pos;

template<typename T, typename T3>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_uniform(T3* gpu_x, T* gpu_x0p, ub32 dims, ub64 seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    T x, y, z;

    randomUniform(&state,&x);
    randomUniform(&state,&y);
    
    if(dims == 3)
    {
        randomUniform(&state,&z);
    }
    

    gpu_x[threadIdx.x].x = fma(x,curr_constStr->box_max[0] - curr_constStr->box_min[0],curr_constStr->box_min[0]);
    gpu_x[threadIdx.x].y = fma(y,curr_constStr->box_max[1] - curr_constStr->box_min[1],curr_constStr->box_min[1]);

    if(dims == 3)
    {
        gpu_x[threadIdx.x].z = fma(z,curr_constStr->box_max[2] - curr_constStr->box_min[2],curr_constStr->box_min[2]);
    }

    gpu_x0p[threadIdx.x] = (T) 1.0f; 
}

template<typename T, typename T3>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_exponential(T3* gpu_x, T* gpu_x0p, ub32 dims, time_t seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    // z is distributed exponentially
    T sigt = curr_constStr->sigt * ((T)2.0);
    T dmax = 0.0;
    if(dims == 3)
    {
        dmax = curr_constStr->box_max[2] - curr_constStr->box_min[2];
    }
    if(dims == 2)
    {
        dmax = curr_constStr->box_max[1] - curr_constStr->box_min[1];
    }

    T rmax = 1-exp(-dmax*sigt);

    T rand_z;
    randomUniform(&state,&rand_z);
    T d = -log(fma(-rand_z,rmax,(T)1.0))/(sigt);

    if(dims == 3)
    {
        gpu_x[threadIdx.x].z = d + curr_constStr->box_min[2];
    }
    if(dims == 2)
    {
        gpu_x[threadIdx.x].y = d + curr_constStr->box_min[1];
    }
    
    gpu_x0p[threadIdx.x] = (exp(-sigt * d) / (rmax)) * sigt * dmax;

    // x and y are distributed uniformly
    T x;
    randomUniform(&state,&x);
    gpu_x[threadIdx.x].x = fma(x,curr_constStr->box_max[0] - curr_constStr->box_min[0],curr_constStr->box_min[0]);

    if(dims == 3)
    {
        T y;
        randomUniform(&state,&y);
        gpu_x[threadIdx.x].y = fma(y,curr_constStr->box_max[1] - curr_constStr->box_min[1],curr_constStr->box_min[1]);
    }
    
}

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub32 sample_position_flag, ub32 is_correlation, ub32 dims, ub64 curr_seed)
{
    if(sample_position_flag == 1)
    {        
        // unifrom sampling
        samplePosition_uniform<T,T3><<<1, THREADS_NUM>>>(x0,px0,dims,curr_seed);
    }
    else if(sample_position_flag == 2)
    {
        samplePosition_exponential<T,T3><<<1, THREADS_NUM>>>(x0,px0,dims,curr_seed);

        // exponential sampling
        if(!is_correlation)
        {
            // in field we take the sqrt of px0
            global_sqrt<T><<<1, THREADS_NUM>>>(px0);
        }
    }
}

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_smp_pos, constMem, sizeof(constStructre<T>));
}

template void samplePosition<double,double3>(double3 *x0, double* px0, ub32 sample_position_flag, ub32 is_correlation, ub32 dims, ub64 curr_seed);
template void samplePosition<float,float3>(float3 *x0, float* px0, ub32 sample_position_flag, ub32 is_correlation, ub32 dims, ub64 curr_seed);

template void samplePosition_setConstMem<double>(constStructre<double> *constMem);
template void samplePosition_setConstMem<float>(constStructre<float> *constMem);
