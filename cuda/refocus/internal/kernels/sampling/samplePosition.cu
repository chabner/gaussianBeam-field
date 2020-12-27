#include "refocus_header.h"
#include "refocus_sampling.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ constStructre<double> constStr_smp_pos;

template<typename T, typename T3>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_uniform(T3* gpu_x, T* px0, ub64 seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    T x, y, z;

    randomUniform(&state,&x);
    randomUniform(&state,&y);
    randomUniform(&state,&z);

    gpu_x[threadIdx.x].x = fma(x,curr_constStr->box_max[0] - curr_constStr->box_min[0],curr_constStr->box_min[0]);
    gpu_x[threadIdx.x].y = fma(y,curr_constStr->box_max[1] - curr_constStr->box_min[1],curr_constStr->box_min[1]);
    gpu_x[threadIdx.x].z = fma(z,curr_constStr->box_max[2] - curr_constStr->box_min[2],curr_constStr->box_min[2]);

    px0[threadIdx.x] = (T) 1.0;
}

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub64 curr_seed)
{
    samplePosition_uniform<T,T3><<<1, THREADS_NUM>>>(x0,px0,curr_seed);
}

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_smp_pos, constMem, sizeof(constStructre<T>));
}

template void samplePosition<double,double3>(double3 *x0, double *px0, ub64 curr_seed);
template void samplePosition<float,float3>(float3 *x0, float *px0, ub64 curr_seed);

template void samplePosition_setConstMem<double>(constStructre<double> *constMem);
template void samplePosition_setConstMem<float>(constStructre<float> *constMem);
