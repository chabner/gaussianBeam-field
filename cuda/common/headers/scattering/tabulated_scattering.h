#pragma once

#include "globalMathKernels.h"

/*
 * The input tabulated_values is a complex function of the amplitude 
 * function (sqrt of phase function), defined from angle [0,2pi] in 2D and 
 * [0,pi] in 3D. tabulated_values[0] is the amplitude function in angle 
 * theta = 0, and tabulated_values[tabulated_entries_num - 1] is related to
 * angle theta = 2pi in 2D and theta = pi in 3D.
 */

template<typename T> struct tabulatedScatter
{
    ub32 tabulated_entries;
    T* evalAmp;
    T* sampleIcdf;
};

template<typename T, typename T2>
void* init_tabulated_scatterer(ub32 tabulated_entries_num, const T2* tabulated_values, bool is3D);

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const tabulatedScatter<T>* scattering_struct, const T3* D, const T3* W);

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const tabulatedScatter<T>* sample_function, ub64 seed, ub32 dims_num);

template<typename T, typename T3, ub32 dims>
__device__ void sample_direction(T3 *w, const tabulatedScatter<T>* sample_function, curandState_t *state);
