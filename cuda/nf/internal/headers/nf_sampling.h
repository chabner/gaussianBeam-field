#pragma once

// I.S. - NOT IMPLEMENTED:
// (1) wavelength dependency (assuming lambda = 1)
// (2) directional beams

#include "nf_header.h"
#include "globalMathKernels.h"
#include "hg_scattering.h"
#include "random_scattering.h"

#define IS_MAX_BLOCKS 131072 // 1024 * 128, corresponds to 1GB of doubles

template<typename T> struct importanceSampling
{
    T* samplePdf; // z0_sample * alpha * samples_num
    T* sampleCdf;
    ub32 z0_sample_num; // number of z0 grid samples
    ub32 z_sample_num; // a multipication of 1024 (i.e. z0_sample_num = 10 => 10240 samples per z0 and alpha)
    ub32* element_to_z0_idx;   // elemNum
    ub32* element_to_z0_idx_2; // elemNum in case of correlation
    T alphaPdf[MAX_DIM];
    T alphaCdf[MAX_DIM];
    T* z0_sample; // z0 values of each sample
    T* z_sample;  // z  values of each sample
    
    // temporal data used by the IS algorithm
    ub32 px_lines; // number of px calculations computed in each iteration before summing to the final probability
    T* px_table; // table stores temporal of px calculations
    ub32* z_idx; // the index in table related to the sampled z
    
};


template<typename T,typename T2,typename T3>
importanceSampling<T>* sampling_preprocess(input_st<T,T2,T3> *data_in, void* gpuGlobalData);

template<typename T>
void sampling_preprocess_free(importanceSampling<T> *IS_struct);

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub32 *n, ub32 *c, ub32 *x_rep,
        ub32 sample_position_flag, ub32 is_correlation, ub32 total_elements,
        importanceSampling<T>* IS_struct, const void* globalMem, ub64 curr_seed, T min_prob);

template<typename T, typename T2, typename T3>
void sampleDirection(const T3 *x0, T3 *w0, T* pw0, const ub32 *n, const ub32 *c,
        ub32 is_correlation, ub32 sample_direction_flag,
        const void* globalMem, ub64 curr_seed, T min_prob);

template<typename T>
void sampling_preprocess_setConstMem(constStructre<T> *constMem, ub32 number_of_elements);

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem, ub32 number_of_elements);

template<typename T>
void sampleDirection_setConstMem(constStructre<T> *constMem, ub32 number_of_elements, bool ipr);
