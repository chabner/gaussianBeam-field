#pragma once

#include "nf_header.h"
#include "sampleScattering.h"
#include "globalMathKernels.h"

template<typename T, typename T2, typename T3>
void singleScattering(T2 *us, const T3 *x0, const T2* constPath, ub32 is_correlation, ub32 total_elements, const void* globalMem);

template<typename T>
void singleScattering_setConstMem(constStructre<T> *constMem);

template<typename T, typename T2, typename T3>
void multipleScattering(T2 *u, const T3 *x0, const T3 *xb, const T3 *w0, const T3 *wb, const T2 *constPath, const pathsList<T>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, const void* globalMem);

template<typename T>
void multipleScattering_setConstMem(constStructre<T> *constMem, ub32 number_of_elements);
