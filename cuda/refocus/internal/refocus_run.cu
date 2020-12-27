#include "refocus_interface.h"
#include "refocus_header.h"
#include "refocus_sampling.h"
#include "refocus_scattering.h"
#include "standardC.h"
#include "refocus_build.h"

#include "sampleScattering.h"

template<typename T,typename T3>
void copy_input_data(T3 *dest, input_data_triplet<T> *in_source, ub32 *dimsArray)
{
    dest->x = in_source->data_ptr_x[c_sub2ind(dimsArray,in_source->data_dims_num_x,in_source->data_size_x)];
    dest->y = in_source->data_ptr_y[c_sub2ind(dimsArray,in_source->data_dims_num_y,in_source->data_size_y)];
    dest->z = in_source->data_ptr_z[c_sub2ind(dimsArray,in_source->data_dims_num_z,in_source->data_size_z)];
}

template<typename T>
void copy_input_data(T *dest, input_data<T> *in_source, ub32 *dimsArray)
{
    *dest = in_source->data_ptr[c_sub2ind(dimsArray,in_source->data_dims_num,in_source->data_size)];
}

template<typename T,typename T2,typename T3>
run_return refocus_run(input_st<T,T2,T3> *data_in, T2 *u, T2 *us, double *norm_factor, ub64 *total_iterations)
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
    T k;
    
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
            copy_input_data<T,T3>(&cpuEntry.pixel_1.illuminationP, &data_in->light_point, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.viewP, &data_in->view_point, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.illuminationDir, &data_in->light_direction, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_1.viewDir, &data_in->view_direction, dimsArray);

            copy_input_data<T>(&cpuEntry.pixel_2.k, &data_in->k_2, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.illuminationP, &data_in->light_point_2, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.viewP, &data_in->view_point_2, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.illuminationDir, &data_in->light_direction_2, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel_2.viewDir, &data_in->view_direction_2, dimsArray);

            *currentEntryStructre = cpuEntry;
            currentEntryStructre++;

        }
        k = cpuEntryStructre->pixel_1.k;
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
            copy_input_data<T,T3>(&cpuEntry.pixel.illuminationP, &data_in->light_point, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel.viewP, &data_in->view_point, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel.illuminationDir, &data_in->light_direction, dimsArray);
            copy_input_data<T,T3>(&cpuEntry.pixel.viewDir, &data_in->view_direction, dimsArray);

            *currentEntryStructre = cpuEntry;
            currentEntryStructre++;
        }
        k = cpuEntryStructre->pixel.k;
        cudaMemcpy(gpuGlobalData, cpuEntryStructre, data_in->total_elements * sizeof(entryStructre_field<T,T3>), cudaMemcpyHostToDevice);
        free(cpuEntryStructre);
    }

    // constant memory
    T V = (data_in->box_max.x - data_in->box_min.x) * 
          (data_in->box_max.y - data_in->box_min.y) *
          (data_in->box_max.z - data_in->box_min.z);

    constStructre<T> *constMem;
    constMem = (constStructre<T> *) malloc(sizeof(constStructre<T>));

    constMem->box_min[0] = data_in->box_min.x;
    constMem->box_min[1] = data_in->box_min.y;
    constMem->box_min[2] = data_in->box_min.z;

    constMem->box_max[0] = data_in->box_max.x;
    constMem->box_max[1] = data_in->box_max.y;
    constMem->box_max[2] = data_in->box_max.z;

    constMem->V = V;

    constMem->sigt = data_in->sigt / 2.0;

    constMem->aperture_kappa_l = data_in->aperture_kappa_l;
    constMem->aperture_kappa_v = data_in->aperture_kappa_v;

    singleScattering_setConstMem<T>(constMem);
    refocus_setConstMem<T>(constMem);
    multipleScattering_setConstMem<T>(constMem);
    samplePosition_setConstMem<T>(constMem);

    free(constMem);

    // scattering
    void *scatterig_function;
    switch(data_in->scattering_input.type)
    {
        case 1:
            scatterig_function = 0;
            break;
        case 2:
            scatterig_function = init_tabulated_scatterer<T,T2>(data_in->scattering_input.tabulated_entries,
                    data_in->scattering_input.tabulated_values, true);
            break;
        case 3:
            scatterig_function = init_hg_scatterer<T>(data_in->scattering_input.g);
            break;
    }
        
    // refocus
    ff_refocus_struct<T2,T3>* ff_refocus = init_refocus<T,T2,T3>(&data_in->refocus, data_in->is_correlation, 
        data_in->total_elements, gpuGlobalData);
            
    if(ff_refocus == 0)
    {           
        return ERR_BAD_MEM_ALLOCATION;
    }


    // allocate memory for the path parameters
    ub32 num_of_paths_inside;
    T3 *x0, *xb, *w0, *wb, *w;
    T *pw0, *px0;
    T2 *const_path;
    pathsList<T>* pl;
    pathsList<T> pl_cpu;

    cudaMalloc(&x0, THREADS_NUM * sizeof(T3));
    cudaMalloc(&xb, THREADS_NUM * sizeof(T3));
    cudaMalloc(&w0, THREADS_NUM * sizeof(T3));
    cudaMalloc(&wb, THREADS_NUM * sizeof(T3));
    cudaMalloc(&w , THREADS_NUM * sizeof(T3));
    
    cudaMalloc(&pw0, THREADS_NUM * sizeof(T));
    cudaMalloc(&px0, THREADS_NUM * sizeof(T));
    cudaMalloc(&const_path, THREADS_NUM * sizeof(T2));
    
    cudaMalloc(&pl, sizeof(pathsList<T>));
    
    ub64 curr_seed = 1000000000ULL * (ub64)(time(0));
    ub32 pL;

    // in random sampling, each iteration is only one path
    ub32 total_paths_num = data_in->refocus.sample_random ? 1 : THREADS_NUM;
    ub32 ms_paths_num;
    
    // main loop
    for(ub32 iter = 0; iter < data_in->iterations_num; iter++)
    {
        samplePosition<T,T3>(x0, px0, curr_seed);
        curr_seed++;

        sampleDirection<T,T3>(w0, pw0, data_in->is_correlation, curr_seed);
        curr_seed++;

        const_path_contribution<T,T2>(const_path, px0, data_in->is_correlation, curr_seed);
        curr_seed++;

        for(ub32 pathNum = 0; pathNum < total_paths_num; pathNum++)
        {
            if(data_in->refocus.sample_random)
            {
                randomize_ff_directions<T,T2,T3>(ff_refocus, &data_in->refocus,
                    data_in->is_correlation, data_in->total_elements, gpuGlobalData, curr_seed);
                curr_seed++;
            }

            singleScattering<T,T2,T3>(us_gpu, x0 + pathNum, ff_refocus, const_path + pathNum, data_in->is_correlation, data_in->total_elements,
                k, data_in->scattering_input.type, scatterig_function);
        }


        pL = 0;
        cudaMemcpy(w , w0, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);
        cudaMemcpy(xb, x0, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);

        px0_mult_pw0<T>(px0, pw0);
        initiate_path_list(pl);
        num_of_paths_inside = THREADS_NUM;
        
        while(1)
        {
            pL++;
            if(pL > 1)
            {
                const_path_contribution<T,T2>(const_path, px0, data_in->is_correlation, curr_seed);
                curr_seed++;

                if(data_in->refocus.sample_random)
                {
                    ms_paths_num = 1;
                }
                else
                {
                    ms_paths_num = pl_cpu.path_count;
                }
                
                for(ub32 pathNum = 0; pathNum < ms_paths_num; pathNum++)
                {
                    multipleScattering<T,T2,T3>(u_gpu, ff_refocus, k,
                        data_in->scattering_input.type, scatterig_function,
                        x0 + pl_cpu.paths_numbers_list[pathNum],
                        xb + pl_cpu.paths_numbers_list[pathNum],
                        w0 + pl_cpu.paths_numbers_list[pathNum],
                        wb + pl_cpu.paths_numbers_list[pathNum],
                        const_path + pl_cpu.paths_numbers_list[pathNum],
                        data_in->is_correlation, data_in->total_elements);
                }
            }
            
            propagate_x<T,T3>(xb, pl, w, data_in->sigt, curr_seed, 3);
            curr_seed++;
            
            
            num_of_paths_inside = is_inside_path_list<T,T3>(pl, xb, data_in->box_min, data_in->box_max, 3);
            cudaMemcpy(&pl_cpu, pl, sizeof(pathsList<T>), cudaMemcpyDeviceToHost);

            if(data_in->refocus.sample_random)
            {
                if(num_of_paths_inside == 0 || pl_cpu.paths_numbers_list[0] != 0)
                {
                    break;
                }
            }
            else
            {
                if(num_of_paths_inside == 0)
                {
                    break;
                }
            }

            cudaMemcpy(wb , w, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);
            
            sample_new_directions<T,T3>(w, data_in->scattering_input.type, scatterig_function, curr_seed, 3);
            curr_seed++;
        }

    }

    // copy the data back to the CPU
    cudaMemcpy(u , u_gpu , data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);
    cudaMemcpy(us, us_gpu, data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);

    if(data_in->refocus.sample_random)
    {
        ub64 paths_num = data_in->refocus.random_directions_number;
        *norm_factor = V * data_in->sigt;

        *total_iterations = (paths_num * paths_num * paths_num * paths_num) * data_in->iterations_num;
    }
    else
    {
        *norm_factor = (data_in->refocus.tabulated_dldv) * (data_in->refocus.tabulated_dldv) * // dlx^2
                       (data_in->refocus.tabulated_dldv) * (data_in->refocus.tabulated_dldv) * // dly^2
                       (data_in->refocus.tabulated_dldv) * (data_in->refocus.tabulated_dldv) * // dvx^2
                       (data_in->refocus.tabulated_dldv) * (data_in->refocus.tabulated_dldv) * // dvy^2
                       V * data_in->sigt;

        *total_iterations = THREADS_NUM * data_in->iterations_num;
    }

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
        
    cudaDeviceReset();
    return VALID_RUN;
}


template run_return refocus_run<double,double2,double3>(input_st<double,double2,double3> *data_in, double2* u, double2* us, double *norm_factor, ub64 *total_iterations);
template run_return refocus_run<float,float2,float3>(input_st<float,float2,float3> *data_in, float2* u, float2* us, double *norm_factor, ub64 *total_iterations);
