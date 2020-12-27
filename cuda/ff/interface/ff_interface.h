#pragma once

#include "standardTypes.h"

#define MAX_INPUT_DIM 32

template<typename T> struct input_data_triplet
{
    const T* data_ptr_x;
    const T* data_ptr_y;
    const T* data_ptr_z;
    
    ub32 data_dims_num_x;
    ub32 data_dims_num_y;
    ub32 data_dims_num_z;
    
    ub32 data_size_x[MAX_INPUT_DIM];
    ub32 data_size_y[MAX_INPUT_DIM];
    ub32 data_size_z[MAX_INPUT_DIM];
};

template<typename T> struct input_data
{
    const T* data_ptr;
    ub32 data_dims_num;
    ub32 data_size[MAX_INPUT_DIM];
};

template<typename T, typename T2> struct input_IS
{
    ub32 position_type;
    ub32 direction_type;
    ub32 is_mean;
    
    // in case of g0 sampling
    T g0;
    
    // in case of tabulated scattering
    ub32 tabulated_entries;
    const T2 *tabulated_values;
};

template<typename T, typename T2> struct input_scattering
{
    // scattering type values:
    // 1 - random
    // 2 - tabulated
    // 3 - HG
    
    ub32 type;
    
    // in case of HG
    T g;
    
    // in case of tabulated
    ub32 tabulated_entries;
    const T2 *tabulated_values;
};

template<typename T, typename T2, typename T3> struct input_st
{
    ub32 cuda_device_num;                                // number of gpu device, from 0 to maximal device num
    ub32 precision;                                      // 0 is single and 1 is double
    ub32 vector_dims_num;                                // 2 for 2D and 3 for 3D
    ub32 is_correlation;                                 // 0 for field and 1 for correlation
    ub32 is_cbs;                                         // is the simulation includes coherent back scattering
    ub32 iterations_num;                                 // number of iterations
    T3 box_min;                                          // minimal size of box
    T3 box_max;                                          // maximal size of box
    T sigt;                                              // 1/MFP
    input_scattering<T,T2> scattering_input;             // scattering type and parameters
    input_IS<T,T2> iS;                                   // importance sampling structre
    input_data<T> k;                                     // wavenumber values of u1
    input_data<T> k_2;                                   // wavenumber values of u2
    input_data_triplet<T> light_direction;               // light point of u1
    input_data_triplet<T> light_attenuation_direction;   // light direction of u1
    input_data_triplet<T> view_direction;                // view point of u1
    input_data_triplet<T> view_attenuation_direction;    // view direction of u1
    input_data_triplet<T> light_direction_2;             // light point of u2
    input_data_triplet<T> light_attenuation_direction_2; // light direction of u2
    input_data_triplet<T> view_direction_2;              // view point of u2
    input_data_triplet<T> view_attenuation_direction_2;  // view direction of u2
    
    // calculated by the algorithm
    ub32 out_dims_num;                                   // number of output dims
    ub32 out_size[MAX_INPUT_DIM];                        // size of output
    ub32 total_elements;                                 // total output elements
};

typedef enum calc_size_return
{
    VALID_SIZE,
    ERR_K,
    ERR_LIGHT_P_X,
    ERR_LIGHT_P_Y,
    ERR_LIGHT_P_Z,
    ERR_LIGHT_D_X,
    ERR_LIGHT_D_Y,
    ERR_LIGHT_D_Z,
    ERR_VIEW_P_X,
    ERR_VIEW_P_Y,
    ERR_VIEW_P_Z,
    ERR_VIEW_D_X,
    ERR_VIEW_D_Y,
    ERR_VIEW_D_Z,
    ERR_LIGHT_P_X_2,
    ERR_LIGHT_P_Y_2,
    ERR_LIGHT_P_Z_2,
    ERR_LIGHT_D_X_2,
    ERR_LIGHT_D_Y_2,
    ERR_LIGHT_D_Z_2,
    ERR_VIEW_P_X_2,
    ERR_VIEW_P_Y_2,
    ERR_VIEW_P_Z_2,
    ERR_VIEW_D_X_2,
    ERR_VIEW_D_Y_2,
    ERR_VIEW_D_Z_2
} calc_size_return;

typedef enum run_return
{
    VALID_RUN,
    ERR_GPU_DEVICE
} run_return;

template<typename T,typename T2,typename T3>
calc_size_return calculate_size(input_st<T,T2,T3> *data_in);

template<typename T,typename T2,typename T3>
run_return ff_run(input_st<T,T2,T3> *data_in, T2* u, T2* us, double *norm_factor, ub64 *total_iterations);
