#pragma once

#include "globalMathKernels.h"

template<typename T> struct HGscatter
{
    T g;
    T g_up;
    T g_down1;
    T g_down2;
};

template<typename T>
void* init_hg_scatterer(T g);

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const HGscatter<T>* scattering_struct, const T3* D, const T3* W);

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const HGscatter<T>* sample_function, ub64 seed, ub32 dims_num);

template<typename T, typename T3, ub32 dims>
__device__ void sample_direction(T3 *w, const HGscatter<T>* sample_function, curandState_t *state);
