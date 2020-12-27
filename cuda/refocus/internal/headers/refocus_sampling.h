#pragma once

#include "refocus_header.h"
#include "globalMathKernels.h"
#include "random_scattering.h"

// only randomized IS is implemented here

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub64 curr_seed);

template<typename T, typename T3>
void sampleDirection(T3 *w0, T* pw0, ub32 is_correlation, ub64 curr_seed);

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem);

