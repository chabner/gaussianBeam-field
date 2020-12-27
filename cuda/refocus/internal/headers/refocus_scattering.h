#pragma once

#include "refocus_header.h"
#include "refocus_build.h"
#include "sampleScattering.h"
#include "globalMathKernels.h"
#include "blasMatrix.h"

template<typename T, typename T2, typename T3>
void singleScattering(T2 *us, const T3 *x0,
	ff_refocus_struct<T2,T3>* ff_refocus,
	const T2* constPath, ub32 is_correlation, ub32 total_elements,
	T k, ub32 scattering_paramenter, const void* scatteringStruct);

template<typename T>
void singleScattering_setConstMem(constStructre<T> *constMem);

template<typename T, typename T2, typename T3>
void multipleScattering(T2 *u, const ff_refocus_struct<T2,T3>* ff_refocus, T k,
                        ub32 scattering_paramenter, const void* scatteringStruct,
                        const T3 *x0, const T3* xb, const T3* w0, const T3* wb, const T2* constPath,
                        ub32 is_correlation, ub32 total_elements);

template<typename T>
void multipleScattering_setConstMem(constStructre<T> *constMem);

