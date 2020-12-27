#pragma once

#include "ff_header.h"
#include "sampleScattering.h"
#include "globalMathKernels.h"
#include "hg_scattering.h"
#include "random_scattering.h"
#include "tabulated_scattering.h"

template<typename T, typename T3>
struct importanceSampling
{
    // in case of g0 hg sampling
    HGscatter<T> *g0_sampling;
    
    // in case of tablated sampling
    tabulatedScatter<T> *f0_sampling;
    
    T3* directions_list;
    ub32 directions_num;
};


template<typename T, typename T2, typename T3>
importanceSampling<T,T3>* sampling_preprocess(const input_IS<T,T2> *iS, const void *gpuGlobalData,
        ub32 is_correlation, bool is3D, ub32 elements_number);

template<typename T, typename T3>
void sampling_preprocess_free(importanceSampling<T,T3> *IS_struct);

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub32 sample_position_flag, ub32 is_correlation, ub32 dims, ub64 curr_seed);

template<typename T, typename T3>
void sampleDirection(T3 *w0, T* pw0, const importanceSampling<T,T3>* iS, ub32 is_correlation,
        ub32 sample_direction_flag, ub32 dims, const void* globalMem, ub64 curr_seed);

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem);
