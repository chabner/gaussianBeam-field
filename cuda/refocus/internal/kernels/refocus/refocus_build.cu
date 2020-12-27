#include "refocus_build.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ constStructre<double> constStr_refocus;

template<typename T, typename T3>
void assign_values(const refocus_data *refocus_data_struct, ff_direction<T3>* cpu_dir_list,
        ub32 current_size, bool is_valid, double x_value, double y_value, double z_value)
{
    cpu_dir_list[current_size].dir.x = (T) x_value;
    cpu_dir_list[current_size].dir.y = (T) y_value;
    cpu_dir_list[current_size].dir.z = (T) z_value;

    if(refocus_data_struct->bias_attenuation)
    {
        cpu_dir_list[current_size].dir_attenuation.x = (T) 0.0;
        cpu_dir_list[current_size].dir_attenuation.y = (T) 0.0;
        cpu_dir_list[current_size].dir_attenuation.z = (T) 1.0;
    }
    else
    {
        cpu_dir_list[current_size].dir_attenuation.x = (T) x_value;
        cpu_dir_list[current_size].dir_attenuation.y = (T) y_value;
        cpu_dir_list[current_size].dir_attenuation.z = (T) z_value;
    }

    cpu_dir_list[current_size].is_active = (ub32) is_valid;
}

template<typename T, typename T3, ub32 correlation_number>
__device__ pixel_entry<T,T3> get_pixel_data(const entryStructre_correlation<T,T3>* globalMem, ub32 pixel_num)
{
    if(correlation_number == 1)
    {
        return (globalMem + pixel_num)->pixel_1;
    }
    else
    {
        return (globalMem + pixel_num)->pixel_2;
    }
}

template<typename T, typename T3, ub32 correlation_number>
__device__ pixel_entry<T,T3> get_pixel_data(const entryStructre_field<T,T3>* globalMem, ub32 pixel_num)
{
    return (globalMem + pixel_num)->pixel;
}

template<typename T, typename T2, typename T3, typename scatteringType, ub32 correlation_number>
__global__ void compute_wv_kernel(ub32 is_binary_aperture, T2* wv, const scatteringType* globalMem, ff_direction<T3> *directions_list, 
        ub32 total_elements, ub32 total_directions)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_refocus;
    ub32 w_elem_num = blockIdx.x * blockDim.x + threadIdx.x;

    if(w_elem_num < (total_elements * total_directions))
    {
        ub32 directional_element = w_elem_num % total_directions;
        ub32 nf_element = w_elem_num / total_directions;

        ff_direction<T3> directional_data = directions_list[directional_element];
        pixel_entry<T,T3> current_pixel = get_pixel_data<T,T3,correlation_number>(globalMem, nf_element);

        if(directional_data.is_active)
        {
            T2 exp_val;
            if(is_binary_aperture)
            {
                exp_val.x = (T) 0.0;
                exp_val.y = (T) current_pixel.k * 
                        (directional_data.dir.x * current_pixel.viewP.x + 
                         directional_data.dir.y * current_pixel.viewP.y +
                         directional_data.dir.z * current_pixel.viewP.z);
                wv[w_elem_num] = complexExponent(exp_val) / abs(directional_data.dir.z);
            }
            else
            {
                exp_val.x = -1 * curr_constStr->aperture_kappa_v + curr_constStr->aperture_kappa_v *
                        (directional_data.dir.x * current_pixel.viewDir.x + 
                         directional_data.dir.y * current_pixel.viewDir.y +
                         directional_data.dir.z * current_pixel.viewDir.z);

                exp_val.y = current_pixel.k * 
                        (directional_data.dir.x * current_pixel.viewP.x + 
                         directional_data.dir.y * current_pixel.viewP.y +
                         directional_data.dir.z * current_pixel.viewP.z);

                if(abs(exp_val.x) < 10)
                {
                    wv[w_elem_num] = ((curr_constStr->aperture_kappa_v / (2 * M_PI)) * complexExponent(exp_val)) / abs(directional_data.dir.z);
                }
                else
                {
                    wv[w_elem_num].x = (T) 0.0;
                    wv[w_elem_num].y = (T) 0.0;
                }
                
            }
        }
        else
        {
            wv[w_elem_num].x = (T) 0.0;
            wv[w_elem_num].y = (T) 0.0;
        }
    }
}

template<typename T, typename T2, typename T3>
void compute_wv(ub32 is_binary_aperture, T2* wv, const void* globalMem, ff_direction<T3> *directions_list, 
        ub32 is_correlation, ub32 corr_num, ub32 total_elements, ub32 total_directions)
{
    ub32 blocks_num = (total_elements * total_directions - 1) / THREADS_NUM + 1;

    if(is_correlation)
    {
        if(corr_num == 1)
        {      
            compute_wv_kernel<T,T2,T3,entryStructre_correlation<T,T3>,1><<<blocks_num,THREADS_NUM>>>
                    (is_binary_aperture, wv, (const entryStructre_correlation<T,T3>*) globalMem, directions_list,
                     total_elements,total_directions);
        }
        if(corr_num == 2)
        {
            compute_wv_kernel<T,T2,T3,entryStructre_correlation<T,T3>,2><<<blocks_num,THREADS_NUM>>>
                    (is_binary_aperture, wv, (const entryStructre_correlation<T,T3>*) globalMem, directions_list,
                     total_elements,total_directions);
        }
    }
    else
    {
        compute_wv_kernel<T,T2,T3,entryStructre_field<T,T3>,1><<<blocks_num,THREADS_NUM>>>
            (is_binary_aperture, wv, (const entryStructre_field<T,T3>*) globalMem, directions_list,
             total_elements,total_directions);
    }
}

template<typename T, typename T2, typename T3, typename scatteringType, ub32 correlation_number>
__global__ void compute_wl_kernel(ub32 is_binary_aperture, T2* wl, const scatteringType* globalMem, ff_direction<T3> *directions_list, 
        ub32 total_elements, ub32 total_directions)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_refocus;
    ub32 w_elem_num = blockIdx.x * blockDim.x + threadIdx.x;

    if(w_elem_num < (total_elements * total_directions))
    {        
        ub32 directional_element = w_elem_num % total_directions;
        ub32 nf_element = w_elem_num / total_directions;

        ff_direction<T3> directional_data = directions_list[directional_element];
        pixel_entry<T,T3> current_pixel = get_pixel_data<T,T3,correlation_number>(globalMem, nf_element);

        if(directional_data.is_active)
        {
            T2 exp_val;
            if(is_binary_aperture)
            {
                exp_val.x = (T) 0.0;

                exp_val.y = (T) -1 * current_pixel.k * 
                        (directional_data.dir.x * current_pixel.illuminationP.x + 
                         directional_data.dir.y * current_pixel.illuminationP.y +
                         directional_data.dir.z * current_pixel.illuminationP.z);
                wl[w_elem_num] = complexExponent(exp_val) / abs(directional_data.dir.z);
            }
            else
            {
                exp_val.x = -1 * curr_constStr->aperture_kappa_l + curr_constStr->aperture_kappa_l *
                        (directional_data.dir.x * current_pixel.illuminationDir.x + 
                         directional_data.dir.y * current_pixel.illuminationDir.y +
                         directional_data.dir.z * current_pixel.illuminationDir.z);
                exp_val.y = -1 * current_pixel.k * 
                        (directional_data.dir.x * current_pixel.illuminationP.x + 
                         directional_data.dir.y * current_pixel.illuminationP.y +
                         directional_data.dir.z * current_pixel.illuminationP.z);

                if(abs(exp_val.x) < 10)
                {
                    wl[w_elem_num] = ((curr_constStr->aperture_kappa_l / (2 * M_PI)) * complexExponent(exp_val)) / abs(directional_data.dir.z);    
                }
                else
                {
                    wl[w_elem_num].x = (T) 0.0;
                    wl[w_elem_num].y = (T) 0.0;
                }

            }
        }
        else
        {
            wl[w_elem_num].x = (T) 0.0;
            wl[w_elem_num].y = (T) 0.0;
        }

    }
}

template<typename T, typename T2, typename T3>
void compute_wl(ub32 is_binary_aperture, T2* wl, const void* globalMem, ff_direction<T3> *directions_list, 
        ub32 is_correlation, ub32 corr_num, ub32 total_elements, ub32 total_directions)
{
    ub32 blocks_num = (total_elements * total_directions - 1) / THREADS_NUM + 1;

    if(is_correlation)
    {
        if(corr_num == 1)
        {
            compute_wl_kernel<T,T2,T3,entryStructre_correlation<T,T3>,1><<<blocks_num,THREADS_NUM>>>
                    (is_binary_aperture, wl, (const entryStructre_correlation<T,T3>*) globalMem, directions_list,
                     total_elements,total_directions);
        }
        if(corr_num == 2)
        {
            compute_wl_kernel<T,T2,T3,entryStructre_correlation<T,T3>,2><<<blocks_num,THREADS_NUM>>>
                    (is_binary_aperture, wl, (const entryStructre_correlation<T,T3>*) globalMem, directions_list,
                     total_elements,total_directions);
        }
    }
    else
    {
        compute_wl_kernel<T,T2,T3,entryStructre_field<T,T3>,1><<<blocks_num,THREADS_NUM>>>
            (is_binary_aperture, wl, (const entryStructre_field<T,T3>*) globalMem, directions_list,
             total_elements,total_directions);
    }
}

template<typename T, typename T2, typename T3>
ff_refocus_struct<T2,T3>* init_refocus(const refocus_data *refocus_data_struct, ub32 is_correlation, 
        ub32 total_elements, const void* globalMem)
{
    ff_refocus_struct<T2,T3>* ff_refocus = 
            (ff_refocus_struct<T2,T3>*) malloc(sizeof(ff_refocus_struct<T2,T3>));

    // init cublas handler
    cublasCreate(&(ff_refocus->cublas_handle));
    
    cublasCreate(&(ff_refocus->cublas_handle_device_alpha));
    cublasSetPointerMode(ff_refocus->cublas_handle_device_alpha, CUBLAS_POINTER_MODE_DEVICE);

    ub32 directions_num = 0;

    // init all directions
    if(refocus_data_struct->sample_random)
    {
        // in random direction sampling, we cannot pre-calculate w
        
        // make sure that the number of sampled directions is divided by 8
        // for benefit tensors cores
        directions_num = refocus_data_struct->random_directions_number;
        if(directions_num % 8 != 0)
        {
            directions_num = directions_num - directions_num % 8 + 8;
        }

        ff_refocus->num_directions = directions_num;
        if(cudaMalloc(&ff_refocus->d_l, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
        {
            return 0;
        }
        if(cudaMalloc(&ff_refocus->d_v, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
        {
            return 0;
        }

        if(is_correlation)
        {
            if(cudaMalloc(&ff_refocus->d_l_2, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
            {
                return 0;
            }
            if(cudaMalloc(&ff_refocus->d_v_2, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
            {
                return 0;
            }
        }
    }
    else
    {
        // we assume double to avoid numerical issues when compute z_value
        double x_value, y_value, z_value;
        ub32 maximal_possible_size = 0;
        ub32 x_y_max_num = ((ub32) 2 * ceil((refocus_data_struct->max_xy_value / refocus_data_struct->tabulated_dldv)) + 1 );

        if(refocus_data_struct->sample_forward)
        {
            maximal_possible_size++;
        }
        if(refocus_data_struct->sample_backward)
        {
            maximal_possible_size++;
        }
        maximal_possible_size *= x_y_max_num * x_y_max_num;
        maximal_possible_size += 8;

        ff_direction<T3>* cpu_dir_list = (ff_direction<T3>*) malloc(maximal_possible_size * sizeof(ff_direction<T3>));

        for(x_value = -1 * refocus_data_struct->max_xy_value; 
            x_value <= refocus_data_struct->max_xy_value;
            x_value += refocus_data_struct->tabulated_dldv)
        {
            for(y_value = -1 * refocus_data_struct->max_xy_value; 
                y_value <= refocus_data_struct->max_xy_value;
                y_value += refocus_data_struct->tabulated_dldv)
            {
                z_value = 1.0 - x_value * x_value - y_value * y_value;

                if(z_value < 0.0 || z_value > 1.0 || x_value < -1.0 || x_value > 1.0 || y_value < -1.0 || y_value > 1.0)
                {
                    continue;
                }

               z_value = sqrt(z_value);

                if(z_value < 0.01)
                {
                    continue;
                }

                // for non zero values, we also have may have negative z values

                if(refocus_data_struct->sample_forward)
                {
                   assign_values<T,T3>(refocus_data_struct,cpu_dir_list,directions_num,true, x_value, y_value, z_value);
                   directions_num++;
                }

                if(refocus_data_struct->sample_backward)
                {
                   assign_values<T,T3>(refocus_data_struct,cpu_dir_list,directions_num,true, x_value, y_value, -1.0 * z_value);
                   directions_num++;
                }

            }
        }

        while(directions_num % 8 != 0)
        {
           assign_values<T,T3>(refocus_data_struct,cpu_dir_list,directions_num,false, 0.0, 0.0, 0.0);
           directions_num++;
        }

        // now we know the number of ff directions
        ff_refocus->num_directions = directions_num;

        if(cudaMalloc(&ff_refocus->d_l, directions_num * sizeof(ff_direction<T3>))!= cudaSuccess)
        {
            return 0;
        }

        if(cudaMalloc(&ff_refocus->d_v, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
        {
            return 0;
        }

        if(cudaMemcpy(ff_refocus->d_l, cpu_dir_list, directions_num * sizeof(ff_direction<T3>), cudaMemcpyHostToDevice) != cudaSuccess)
        {
            return 0;
        }

        if(cudaMemcpy(ff_refocus->d_v, ff_refocus->d_l, directions_num * sizeof(ff_direction<T3>), cudaMemcpyDeviceToDevice) != cudaSuccess)
        {
            return 0;
        }

        if(is_correlation)
        {
            if(cudaMalloc(&ff_refocus->d_l_2, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
            {
                return 0;
            }

            if(cudaMalloc(&ff_refocus->d_v_2, directions_num * sizeof(ff_direction<T3>)) != cudaSuccess)
            {
                return 0;
            }

            if(cudaMemcpy(ff_refocus->d_l_2, ff_refocus->d_l, directions_num * sizeof(ff_direction<T3>), cudaMemcpyDeviceToDevice) != cudaSuccess)
            {
                return 0;
            }

            if(cudaMemcpy(ff_refocus->d_v_2, ff_refocus->d_l, directions_num * sizeof(ff_direction<T3>), cudaMemcpyDeviceToDevice) != cudaSuccess)
            {
                return 0;
            }
        }
                
        free(cpu_dir_list);

    }

    if(cudaMalloc(&ff_refocus->w_l, directions_num * total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->w_v, directions_num * total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(is_correlation)
    {
        if(cudaMalloc(&ff_refocus->w_l_2, directions_num * total_elements * sizeof(T2)) != cudaSuccess)
        {
            return 0;
        }

        if(cudaMalloc(&ff_refocus->w_v_2, directions_num * total_elements * sizeof(T2)) != cudaSuccess)
        {
            return 0;
        }

        if(cudaMalloc(&ff_refocus->mult_res_2, total_elements * sizeof(T2)) != cudaSuccess)
        {
            return 0;
        }
    }

    if(cudaMalloc(&ff_refocus->single_scattering_ff, directions_num * directions_num * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->multiple_scattering_view_ff, directions_num * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->multiple_scattering_illumination_ff, directions_num * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->mult_w_v_single_scattering, directions_num * total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->mult_w_v_multiple_scattering, total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->mult_w_l_multiple_scattering, total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }
    
    if(cudaMalloc(&ff_refocus->mult_res, total_elements * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    if(cudaMalloc(&ff_refocus->ones_vector, directions_num * sizeof(T2)) != cudaSuccess)
    {
        return 0;
    }

    T2 *cpu_ones_dir_list = (T2*) malloc(sizeof(T2) * directions_num);
    for(ub32 iter_num = 0; iter_num < directions_num; iter_num++)
    {
        cpu_ones_dir_list[iter_num].x = (T) 1.0;
        cpu_ones_dir_list[iter_num].y = (T) 0.0;
    }

    if(cudaMemcpy(ff_refocus->ones_vector, cpu_ones_dir_list, directions_num * sizeof(T2), cudaMemcpyHostToDevice) != cudaSuccess)
    {
        return 0;
    }

    free(cpu_ones_dir_list);

    // for tabulated MC, we can compute now the w weights
    if(!refocus_data_struct->sample_random)
    {        
        compute_wv<T,T2,T3>(refocus_data_struct->binary_aperture, ff_refocus->w_v, globalMem, ff_refocus->d_v, 
            is_correlation, 1, total_elements, directions_num);
        compute_wl<T,T2,T3>(refocus_data_struct->binary_aperture, ff_refocus->w_l, globalMem, ff_refocus->d_l, 
            is_correlation, 1, total_elements, directions_num);

            if(is_correlation)
            {
                compute_wv<T,T2,T3>(refocus_data_struct->binary_aperture, ff_refocus->w_v_2, globalMem, ff_refocus->d_v_2, 
                    is_correlation, 2, total_elements, directions_num);
                compute_wl<T,T2,T3>(refocus_data_struct->binary_aperture, ff_refocus->w_l_2, globalMem, ff_refocus->d_l_2, 
                    is_correlation, 2, total_elements, directions_num);
            }
    }

    return ff_refocus;
}

template<typename T, typename T2, typename T3>
__global__ void randomize_ff_directions_kernel(ub32 directions_num, ff_direction<T3> *directions_list,
        ub32 bias_attenuation, ub64 seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_refocus;
    ub32 w_elem_num = blockIdx.x * blockDim.x + threadIdx.x;

    curandState_t state;
    curand_init(seed, w_elem_num, 0, &state);

    if(w_elem_num < directions_num)
    {
        T3 center;
        center.x = (T) 0.0;
        center.y = (T) 0.0;
        center.z = (T) 1.0;

        directions_list[w_elem_num].dir = random_vMF_direction<T,T2,T3>(center, curr_constStr->aperture_kappa_l, &state);
        if(bias_attenuation)
        {
            directions_list[w_elem_num].dir_attenuation = center;
        }
        else
        {
            directions_list[w_elem_num].dir_attenuation = directions_list[w_elem_num].dir;
        }
        directions_list[w_elem_num].is_active = true;
    }
}

template<typename T, typename T2, typename T3>
void randomize_ff_directions(ff_refocus_struct<T2,T3> *ff_refocus, const refocus_data *refocus_data_struct,
        ub32 is_correlation, ub32 total_elements, const void* globalMem, ub64 seed)
{
    ub32 total_directions = ff_refocus->num_directions;
    ub32 blocks_num = (total_directions - 1) / 64 + 1;

    randomize_ff_directions_kernel<T,T2,T3><<<blocks_num,64>>>(total_directions, ff_refocus->d_l,
        refocus_data_struct->bias_attenuation, seed);

    randomize_ff_directions_kernel<T,T2,T3><<<blocks_num,64>>>(total_directions, ff_refocus->d_v,
        refocus_data_struct->bias_attenuation, seed);

    if(is_correlation)
    {
        randomize_ff_directions_kernel<T,T2,T3><<<blocks_num,64>>>(total_directions, ff_refocus->d_l_2,
            refocus_data_struct->bias_attenuation, seed);

        randomize_ff_directions_kernel<T,T2,T3><<<blocks_num,64>>>(total_directions, ff_refocus->d_v_2,
            refocus_data_struct->bias_attenuation, seed);
    }

    compute_wv<T,T2,T3>(true, ff_refocus->w_v, globalMem, ff_refocus->d_v, 
        is_correlation, 1, total_elements, total_directions);
    compute_wl<T,T2,T3>(true, ff_refocus->w_l, globalMem, ff_refocus->d_l, 
        is_correlation, 1, total_elements, total_directions);

    if(is_correlation)
    {
        compute_wv<T,T2,T3>(true, ff_refocus->w_v_2, globalMem, ff_refocus->d_v_2, 
            is_correlation, 2, total_elements, total_directions);
        compute_wl<T,T2,T3>(true, ff_refocus->w_l_2, globalMem, ff_refocus->d_l_2, 
            is_correlation, 2, total_elements, total_directions);
    }

}

template<typename T>
void refocus_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_refocus, constMem, sizeof(constStructre<T>));
}

template void refocus_setConstMem<double>(constStructre<double> *constMem);
template void refocus_setConstMem<float>(constStructre<float> *constMem);

template ff_refocus_struct<double2,double3>* init_refocus<double,double2,double3>(const refocus_data *refocus_data_struct, ub32 is_correlation, 
        ub32 total_elements, const void* globalMem);
template ff_refocus_struct<float2,float3>* init_refocus<float,float2,float3>(const refocus_data *refocus_data_struct, ub32 is_correlation, 
        ub32 total_elements, const void* globalMem);

template void randomize_ff_directions<double,double2,double3>(ff_refocus_struct<double2,double3> *ff_refocus, const refocus_data *refocus_data_struct,
        ub32 is_correlation, ub32 total_elements, const void* globalMem, ub64 seed);
template void randomize_ff_directions<float,float2,float3>(ff_refocus_struct<float2,float3> *ff_refocus, const refocus_data *refocus_data_struct,
        ub32 is_correlation, ub32 total_elements, const void* globalMem, ub64 seed);
        