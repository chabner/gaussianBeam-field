#pragma once

#include "hg_scattering.h"
#include "random_scattering.h"
#include "tabulated_scattering.h"
#include "standardTypes.h"
#include "globalMathKernels.cu"

template<typename T>
struct pathsList
{
	ub32 paths_inside[32];
    ub32 paths_numbers_list[THREADS_NUM];
    ub32 path_count;
    T path_length[THREADS_NUM];
};

template<typename T, typename T2>
void const_path_contribution(T2 *const_path, const T* px0, ub32 is_correlation, ub64 curr_seed);

template<typename T, typename T3>
void propagate_x(T3 *x, pathsList<T> *pl, const T3* w, T sigt, ub64 seed, ub32 dims_num);

template<typename T>
void initiate_path_list(pathsList<T> *pl);

template<typename T, typename T3>
ub32 is_inside_path_list(pathsList<T> *pl, const T3 *x, T3 box_min, T3 box_max, ub32 dims_num);

template<typename T, typename T3>
void sample_new_directions(T3 *w, ub32 sample_flag, const void *sample_function, ub64 curr_seed, ub32 dims_num );

template<typename T>
void px0_mult_pw0(T *px0, const T* pw0);
