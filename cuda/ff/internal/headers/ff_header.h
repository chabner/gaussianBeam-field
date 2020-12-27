#pragma once

#define MAX_DIM 32

#include "cuda_runtime.h"
#include "curand.h"
#include "curand_kernel.h"
#include "device_launch_parameters.h"
#include "math.h"
#include "time.h"
#include "standardTypes.h"
#include "ff_interface.h"

template<typename T, typename T3> struct pixel_entry
{
    T3 illumination_direction;
    T3 illumination_attenuation_direction;
    T3 view_direction;
    T3 view_attenuation_direction;
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
};

