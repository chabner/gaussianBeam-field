#pragma once

#include "standardTypes.h"
#include "curand_kernel.h"
#include "curand.h"

template<typename T>
__global__ void global_sqrt(T* x);

template<typename T>
__global__ void global_mult(T* x, const T* y);

template<typename T>
__global__ void global_add(T* x, const T* y);

template<typename T>
__global__ void global_min(T* x, T min_val);

template<typename T2>
__global__ void elementwise_conj_mult(T2* x, const T2* y, ub32 n);

template<typename T2>
__global__ void elementwise_mult(T2* x, const T2* y, ub32 n);

template<typename T, typename T3, ub32 dims>
__device__ void rotate_by_theta(T3 *w, T costheta, curandState_t *state);

template<typename T,typename T3, ub32 dims>
__global__ void normalize_dir(T3 *w);

template<typename T>
__device__ ub32 binary_search(T *search_array, T search_target, ub32 elements_num);
