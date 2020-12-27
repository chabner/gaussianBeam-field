#pragma once

#include "globalMathKernels.h"

typedef ub32 randomScatter;

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const randomScatter* scattering_struct, const T3* D, const T3* W);

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const randomScatter* sample_function, ub64 seed, ub32 dims_num);

template<typename T, typename T3, ub32 dims>
__device__ void sample_direction(T3 *w, const randomScatter* sample_function, curandState_t *state);
