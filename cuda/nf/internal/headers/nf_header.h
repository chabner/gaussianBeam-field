#pragma once

#define MAX_DIM 13

#include "cuda_runtime.h"
#include "curand.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "time.h"
#include "standardTypes.h"
#include "nf_interface.h"

template<typename T, typename T3> struct pixel_entry
{
    T3 illuminationP;
    T3 illuminationDir;
    T3 viewP;
    T3 viewDir;
    T k;
};

template<typename T, typename T3> struct entryStructre_field
{
    pixel_entry<T,T3> pixel;
};

template<typename T, typename T3> struct entryStructre_correlation
{
    pixel_entry<T,T3> pixel_1;
    pixel_entry<T,T3> pixel_2;
};

template<typename T> struct constStructre
{
    T box_min[3];
    T box_max[3];
    T V;
    T sigt; // sigt/2
	T mixtureMu[MAX_DIM];
	T mixtureC[MAX_DIM];
	ub32 mixturesNum;
	T aperture_kappa_l;
	T aperture_kappa_v;
	T aperture_C_l_plus_aperture_C_v_plus_LOG_2_PI;
};

