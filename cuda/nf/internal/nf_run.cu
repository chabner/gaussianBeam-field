#include "nf_interface.h"
#include "nf_header.h"
#include "nf_sampling.h"
#include "nf_scattering.h"
#include "standardC.h"

#include "sampleScattering.h"
#include "kth_smaller.cpp"

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
run_return nf_run(input_st<T,T2,T3> *data_in, T2 *u, T2 *us, double *norm_factor, ub64 *total_iterations, T* min_px0, T* min_pw0)
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
        cudaMemcpy(gpuGlobalData, cpuEntryStructre, data_in->total_elements * sizeof(entryStructre_field<T,T3>), cudaMemcpyHostToDevice);
        free(cpuEntryStructre);
    }

    // determinate if sampling direction from the same beam which sampled for first position
    bool use_same_beam = (data_in->iS.is_same_beam);

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
    constMem->mixturesNum = data_in->scattering_vMF_mixture.mixture_size;

    for(ub32 mixNum = 0; mixNum < constMem->mixturesNum; mixNum++)
    {
        constMem->mixtureMu[mixNum] = data_in->scattering_vMF_mixture.mixture_mu[mixNum];
        constMem->mixtureC[mixNum] = data_in->scattering_vMF_mixture.mixture_c[mixNum] + log(data_in->scattering_vMF_mixture.mixture_alpha[mixNum]);
    }

    constMem->aperture_kappa_l = data_in->aperture_kappa_l;
    constMem->aperture_kappa_v = data_in->aperture_kappa_v;
    constMem->aperture_C_l_plus_aperture_C_v_plus_LOG_2_PI = data_in->aperture_c_l + data_in->aperture_c_v + 2*log(2*M_PI);

    singleScattering_setConstMem<T>(constMem);
    multipleScattering_setConstMem<T>(constMem, data_in->total_elements);
    sampling_preprocess_setConstMem<T>(constMem, data_in->total_elements);
    samplePosition_setConstMem<T>(constMem, data_in->total_elements);
    sampleDirection_setConstMem<T>(constMem, data_in->total_elements, use_same_beam);

    free(constMem);

    // sampler
    importanceSampling<T>* IS_struct;
    IS_struct = sampling_preprocess<T,T2,T3>(data_in,gpuGlobalData);

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
    
    // allocate memory for the path parameters
    ub32 num_of_paths_inside;
    T3 *x0, *xb, *w0, *wb, *w;
    T *px0, *pw0;
    T2 *const_path;
    ub32 *n, *c, *x_rep;
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
    
    cudaMalloc(&n, THREADS_NUM * sizeof(ub32));
    cudaMalloc(&c, THREADS_NUM * sizeof(ub32));
    cudaMalloc(&x_rep, sizeof(ub32));
    
    cudaMemset(x_rep, 0, sizeof(ub32));
    
    ub64 curr_seed = 1000000000ULL * (ub64)(time(0));
    ub32 pL;

    //calculate minimal probability
    if(!data_in->iS.is_min_known)
    {
        ub32 min_k = (ub32) round((data_in->iS.sample_rounds * THREADS_NUM) * data_in->iS.min_percent);
        T *px0_tmp_gpu, *pw0_tmp_gpu;
        cudaMalloc(&px0_tmp_gpu, data_in->iS.sample_rounds * THREADS_NUM * sizeof(T));
        cudaMalloc(&pw0_tmp_gpu, data_in->iS.sample_rounds * THREADS_NUM * sizeof(T));
        ub32 array_pointer = 0;

        for(ub32 round_num = 0; round_num < data_in->iS.sample_rounds; round_num++)
        {
            samplePosition<T,T3>(x0, px0_tmp_gpu + array_pointer, n, c, x_rep,
                data_in->iS.position_type, data_in->is_correlation, data_in->total_elements,
                IS_struct, gpuGlobalData, curr_seed, (T) 0.0);
            curr_seed++;

            sampleDirection<T,T2,T3>(x0, w0, pw0_tmp_gpu + array_pointer, n, c,
                data_in->is_correlation, data_in->iS.direction_type, gpuGlobalData, curr_seed, (T) 0.0);
            curr_seed++;

            array_pointer += THREADS_NUM;
        }

        T *px0_tmp_cpu, *pw0_tmp_cpu;
        px0_tmp_cpu = (T*) malloc(data_in->iS.sample_rounds * THREADS_NUM * sizeof(T));
        pw0_tmp_cpu = (T*) malloc(data_in->iS.sample_rounds * THREADS_NUM * sizeof(T));

        cudaMemcpy(px0_tmp_cpu , px0_tmp_gpu, data_in->iS.sample_rounds * THREADS_NUM * sizeof(T), cudaMemcpyDeviceToHost);
        cudaMemcpy(pw0_tmp_cpu , pw0_tmp_gpu, data_in->iS.sample_rounds * THREADS_NUM * sizeof(T), cudaMemcpyDeviceToHost);

        data_in->iS.min_px0 = kthSmallest<T>(px0_tmp_cpu, data_in->iS.sample_rounds * THREADS_NUM, min_k);
        data_in->iS.min_pw0 = kthSmallest<T>(pw0_tmp_cpu, data_in->iS.sample_rounds * THREADS_NUM, min_k);

        cudaFree(px0_tmp_gpu);
        cudaFree(pw0_tmp_gpu);
        free(px0_tmp_cpu);
        free(pw0_tmp_cpu);
    }

    // main loop
    for(ub32 iter = 0; iter < data_in->iterations_num; iter++)
    {
        samplePosition<T,T3>(x0, px0, n, c, x_rep,
            data_in->iS.position_type, data_in->is_correlation, data_in->total_elements,
            IS_struct, gpuGlobalData, curr_seed, data_in->iS.min_px0);
        curr_seed++;

        sampleDirection<T,T2,T3>(x0, w0, pw0, n, c, data_in->is_correlation,
            data_in->iS.direction_type, gpuGlobalData, curr_seed, data_in->iS.min_pw0);
        curr_seed++;

        const_path_contribution<T,T2>(const_path, px0, data_in->is_correlation, curr_seed);
        curr_seed++;

        singleScattering<T,T2,T3>(us_gpu, x0, const_path, data_in->is_correlation, data_in->total_elements, gpuGlobalData);

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
                
                multipleScattering<T,T2,T3>(u_gpu, x0, xb, w0, wb, const_path, pl, 
                    data_in->is_correlation, data_in->total_elements, num_of_paths_inside, gpuGlobalData);
            }
            
            propagate_x<T,T3>(xb, pl, w, data_in->sigt, curr_seed, 3);
            curr_seed++;
            
            
            if((num_of_paths_inside = is_inside_path_list<T,T3>(pl, xb, data_in->box_min, data_in->box_max, 3)) == 0 )
            {
                break;
            }

            cudaMemcpy(wb , w, THREADS_NUM * sizeof(T3), cudaMemcpyDeviceToDevice);
            
            sample_new_directions<T,T3>(w, data_in->scattering_input.type, scatterig_function, curr_seed, 3);
            curr_seed++;
        }
    }

    // copy the data back to the CPU
    cudaMemcpy(u , u_gpu , data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);
    cudaMemcpy(us, us_gpu, data_in->total_elements * sizeof(T2), cudaMemcpyDeviceToHost);
    *norm_factor = V * data_in->sigt;

    ub32 x_rep_cpu;
    cudaMemcpy(&x_rep_cpu , x_rep , sizeof(ub32), cudaMemcpyDeviceToHost);
    *total_iterations = THREADS_NUM * data_in->iterations_num + x_rep_cpu;

    *min_px0 = data_in->iS.min_px0;
    *min_pw0 = data_in->iS.min_pw0;

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
//     cudaFree(n);
//     cudaFree(c);
//     cudaFree(x_rep);
//     cudaFree(scatterig_function);
// 
//     isampling_preprocess_free<T>(IS_struct);;

    cudaDeviceReset();
    return VALID_RUN;
}

template run_return nf_run<double,double2,double3>(input_st<double,double2,double3> *data_in, double2* u, double2* us, double *norm_factor, ub64 *total_iterations, double *min_px0, double *min_pw0);
template run_return nf_run<float,float2,float3>(input_st<float,float2,float3> *data_in, float2* u, float2* us, double *norm_factor, ub64 *total_iterations, float *min_px0, float *min_pw0);
