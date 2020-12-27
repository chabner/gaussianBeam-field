#include "ff_interface.h"
#include "ff_header.h"
#include "ff_sampling.h"
#include "ff_scattering.h"
#include "standardC.h"

#include "sampleScattering.h"

template<typename T,typename T3>
void copy_input_data(T3 *dest, input_data_triplet<T> *in_source, ub32 *dimsArray, bool is3D)
{
    dest->x = in_source->data_ptr_x[c_sub2ind(dimsArray,in_source->data_dims_num_x,in_source->data_size_x)];
    dest->y = in_source->data_ptr_y[c_sub2ind(dimsArray,in_source->data_dims_num_y,in_source->data_size_y)];
    if(is3D)
    {
        dest->z = in_source->data_ptr_z[c_sub2ind(dimsArray,in_source->data_dims_num_z,in_source->data_size_z)];
    }
}

template<typename T>
void copy_input_data(T *dest, input_data<T> *in_source, ub32 *dimsArray)
{
    *dest = in_source->data_ptr[c_sub2ind(dimsArray,in_source->data_dims_num,in_source->data_size)];
}

template<typename T,typename T2,typename T3>
run_return ff_run(input_st<T,T2,T3> *data_in, T2 *u, T2 *us, double *norm_factor, ub64 *total_iterations)
{
    // set cuda device
    if(cudaSetDevice(data_in->cuda_device_num) != cudaSuccess)
    {
        return ERR_GPU_DEVICE;
    }
    
    // allocate and initiate the output
    T2 *us_gpu, *u_gpu;

    cudaMalloc(&us_gpu, data_in->total_elements * sizeof(T2));
    cudaMemset(us_gpu, 0, data_in->total_elements * sizeof(T2));

    cudaMalloc(&u_gpu, data_in->total_elements * sizeof(T2));
    cudaMemset(u_gpu, 0, data_in->total_elements * sizeof(T2));
    
    // allocate and initiate the global data (all pixel data)
    void *gpuGlobalData;
    ub32 dimsArray[MAX_DIM];
    bool is3D = (data_in->vector_dims_num == 3);
    
    if(data_in->is_correlation)
    {
        cudaMalloc(&gpuGlobalData, data_in->total_elements * sizeof(entryStructre_correlation<T,T3>));

        entryStructre_correlation<T,T3> *cpuEntryStructre, *currentEntryStructre, cpuEntry;
        cpuEntryStructre = (entryStructre_correlation<T,T3> *) malloc(data_in->total_elements * sizeof(entryStructre_correlation<T,T3>));
        currentEntryStructre = cpuEntryStructre;

        for(ub32 elemNum = 0 ; elemNum < data_in->total_elements ; elemNum++ )
        {
            c_ind2sub(dimsArray, data_in->out_size, data_in->out_dims_num, elemNum);

            copy_input_data<T>(&cpuEntry.pixel_1.k, &data_in->k, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.illumination_direction, &data_in->light_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.view_direction, &data_in->view_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.illumination_attenuation_direction, &data_in->light_attenuation_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.view_attenuation_direction, &data_in->view_attenuation_direction, dimsArray, is3D);

            copy_input_data<T>(&cpuEntry.pixel_2.k, &data_in->k_2, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.illumination_direction, &data_in->light_direction_2, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.view_direction, &data_in->view_direction_2, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.illumination_attenuation_direction, &data_in->light_attenuation_direction_2, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.view_attenuation_direction, &data_in->view_attenuation_direction_2, dimsArray, is3D);

            *currentEntryStructre = cpuEntry;
            currentEntryStructre++;
        }
        cudaMemcpy(gpuGlobalData, cpuEntryStructre, data_in->total_elements * sizeof(entryStructre_correlation<T,T3>), cudaMemcpyHostToDevice);
        free(cpuEntryStructre);
    }
    else
    {
        cudaMalloc(&gpuGlobalData, data_in->total_elements * sizeof(entryStructre_field<T,T3>));

        entryStructre_field<T,T3> *cpuEntryStructre, *currentEntryStructre, cpuEntry;
        cpuEntryStructre = (entryStructre_field<T,T3> *) malloc(data_in->total_elements * sizeof(entryStructre_field<T,T3>));
        currentEntryStructre = cpuEntryStructre;

        for(ub32 elemNum = 0 ; elemNum < data_in->total_elements ; elemNum++ )
        {
            c_ind2sub(dimsArray, data_in->out_size, data_in->out_dims_num, elemNum);

            copy_input_data<T>(&cpuEntry.pixel.k, &data_in->k, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel.illumination_direction, &data_in->light_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel.view_direction, &data_in->view_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel.illumination_attenuation_direction, &data_in->light_attenuation_direction, dimsArray, is3D);
            copy_input_data<T,T3>(&cpuEntry.pixel.view_attenuation_direction, &data_in->view_attenuation_direction, dimsArray, is3D);

            *currentEntryStructre = cpuEntry;
            currentEntryStructre++;
        }
        cudaMemcpy(gpuGlobalData, cpuEntryStructre, data_in->total_elements * sizeof(entryStructre_field<T,T3>), cudaMemcpyHostToDevice);
        free(cpuEntryStructre);
    }

    // constant memory
    T V;
    if(is3D)
    {
        V = (data_in->box_max.x - data_in->box_min.x) * 
            (data_in->box_max.y - data_in->box_min.y) *
            (data_in->box_max.z - data_in->box_min.z);
    }
    else
    {
        V = (data_in->box_max.x - data_in->box_min.x) * 
            (data_in->box_max.y - data_in->box_min.y);
    }
    
    constStructre<T> *constMem;
    constMem = (constStructre<T> *) malloc(sizeof(constStructre<T>));

    constMem->box_min[0] = data_in->box_min.x;
    constMem->box_min[1] = data_in->box_min.y;

    constMem->box_max[0] = data_in->box_max.x;
    constMem->box_max[1] = data_in->box_max.y;

    if(is3D)
    {
        constMem->box_min[2] = data_in->box_min.z;
        constMem->box_max[2] = data_in->box_max.z;
    }

    constMem->V = V;

    constMem->sigt = data_in->sigt / 2.0;

    singleScattering_setConstMem<T>(constMem);
    multipleScattering_setConstMem<T>(constMem, data_in->total_elements);
    samplePosition_setConstMem<T>(constMem);

    free(constMem);

    // sampler
    importanceSampling<T,T3>* IS_struct;
    IS_struct = sampling_preprocess<T,T2,T3>(&data_in->iS, gpuGlobalData, data_in->is_correlation, is3D, data_in->total_elements);

    // scattering
    void *scatterig_function;
    switch(data_in->scattering_input.type)
    {
        case 1:
            scatterig_function = 0;
            break;
        case 2:
            scatterig_function = init_tabulated_scatterer<T,T2>(data_in->scattering_input.tabulated_entries,
                    data_in->scattering_input.tabulated_values, is3D);
            break;
        case 3:
            scatterig_function = init_hg_scatterer<T>(data_in->scattering_input.g);
            break;
    }
    
    // allocate memory for the path parameters
    ub32 num_of_paths_inside;
    T3 *x0, *xb, *w0, *wb, *w;
    T *px0, *pw0;
    T2 *const_path;
    pathsList<T>* pl;

    cudaMalloc(&x0, THREADS_NUM * sizeof(T3));
    cudaMalloc(&xb, THREADS_NUM * sizeof(T3));
    cudaMalloc(&w0, THREADS_NUM * sizeof(T3));
    cudaMalloc(&wb, THREADS_NUM * sizeof(T3));
    cudaMalloc(&w , THREADS_NUM * sizeof(T3));
    
    cudaMalloc(&px0, THREADS_NUM * sizeof(T));
    cudaMalloc(&pw0, THREADS_NUM * sizeof(T));
    
    cudaMalloc(&const_path, THREADS_NUM * sizeof(T2));
    
    cudaMalloc(&pl, sizeof(pathsList<T>));
    
    ub64 curr_seed = 1000000000ULL * (ub64)(time(0));
    ub32 pL;
    
    // main loop
    for(ub32 iter = 0; iter < data_in->iterations_num; iter++)
    {
        samplePosition<T,T3>(x0, px0, data_in->iS.position_type, data_in->is_correlation, data_in->vector_dims_num, curr_seed);
        curr_seed++;

        sampleDirection<T,T3>(w0, pw0, IS_struct, data_in->is_correlation, data_in->iS.direction_type, data_in->vector_dims_num, 
                gpuGlobalData, curr_seed);
        curr_seed++;

        const_path_contribution<T,T2>(const_path, px0, data_in->is_correlation, curr_seed);
        curr_seed++;

        singleScattering<T,T2,T3>(us_gpu, x0, const_path, data_in->is_correlation, data_in->total_elements, 
                data_in->vector_dims_num, gpuGlobalData, data_in->scattering_input.type, scatterig_function);

        pL = 0;
        cudaMemcpy(w , w0, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);
        cudaMemcpy(xb, x0, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);

        px0_mult_pw0<T>(px0, pw0);
        initiate_path_list<T>(pl);
        num_of_paths_inside = THREADS_NUM;
        
        while(1)
        {
            pL++;
            if(pL > 1)
            {
                const_path_contribution<T,T2>(const_path, px0, data_in->is_correlation, curr_seed);
                curr_seed++;
                
                multipleScattering<T,T2,T3>(u_gpu, x0, xb, w0, wb, const_path, pl, 
                    data_in->is_correlation, data_in->total_elements, num_of_paths_inside, data_in->vector_dims_num, data_in->is_cbs, 
                    gpuGlobalData, data_in->scattering_input.type, scatterig_function);
            }
            
            propagate_x<T,T3>(xb, pl, w, data_in->sigt, curr_seed, data_in->vector_dims_num);
            curr_seed++;
            
            
            if((num_of_paths_inside = is_inside_path_list<T,T3>(pl, xb, data_in->box_min, data_in->box_max, data_in->vector_dims_num)) == 0 )
            {
                break;
            }

            cudaMemcpy(wb , w, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);
            
            sample_new_directions<T,T3>(w, data_in->scattering_input.type, scatterig_function, curr_seed, data_in->vector_dims_num);
            curr_seed++;
        }
    }

    // copy the data back to the CPU
    cudaMemcpy(u , u_gpu , data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);
    cudaMemcpy(us, us_gpu, data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);
    *norm_factor = V * data_in->sigt;

    *total_iterations = THREADS_NUM * data_in->iterations_num;

//     cudaFree(u_gpu);
//     cudaFree(us_gpu);
//     cudaFree(gpuGlobalData);
//     cudaFree(x0);
//     cudaFree(xb);
//     cudaFree(w0);
//     cudaFree(wb);
//     cudaFree(w);
//     cudaFree(px0);
//     cudaFree(pw0);
//     cudaFree(const_path);
//     cudaFree(pl);
//     cudaFree(scatterig_function);
// 
//     sampling_preprocess_free<T,T3>(IS_struct);

    cudaDeviceReset();
    return VALID_RUN;
}

template run_return ff_run<double,double2,double3>(input_st<double,double2,double3> *data_in, double2* u, double2* us, double *norm_factor, ub64 *total_iterations);
template run_return ff_run<float,float2,float3>(input_st<float,float2,float3> *data_in, float2* u, float2* us, double *norm_factor, ub64 *total_iterations);
