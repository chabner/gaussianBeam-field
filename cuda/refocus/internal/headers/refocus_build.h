#pragma once

#include "refocus_header.h"

template<typename T3> struct ff_direction
{
    T3 dir;
    T3 dir_attenuation;
    ub32 is_active;
};

template<typename T2, typename T3> struct ff_refocus_struct
{
    ff_direction<T3> *d_l;
    ff_direction<T3> *d_v;
    ff_direction<T3> *d_l_2;
    ff_direction<T3> *d_v_2;
    ub32 num_directions;
    
    cublasHandle_t cublas_handle;
    cublasHandle_t cublas_handle_device_alpha;
    
    T2 *w_l;
    T2 *w_v;
    T2 *w_l_2;
    T2 *w_v_2;
    
    T2 *single_scattering_ff;
    T2 *multiple_scattering_view_ff;
    T2 *multiple_scattering_illumination_ff;
    
    T2 *mult_w_v_single_scattering;
    
    T2 *mult_w_v_multiple_scattering;
    T2 *mult_w_l_multiple_scattering;
    
    T2 *mult_res;
    T2 *mult_res_2;
    
    T2 *ones_vector;
};

template<typename T>
void refocus_setConstMem(constStructre<T> *constMem);

template<typename T, typename T2, typename T3>
ff_refocus_struct<T2,T3>* init_refocus(const refocus_data *refocus_data_struct, ub32 is_correlation, 
        ub32 total_elements, const void* globalMem);

template<typename T, typename T2, typename T3>
void randomize_ff_directions(ff_refocus_struct<T2,T3> *ff_refocus, const refocus_data *refocus_data_struct,
        ub32 is_correlation, ub32 total_elements, const void* globalMem, ub64 seed);
