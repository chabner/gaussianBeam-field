#include "ff_header.h"
#include "ff_sampling.h"

template<typename T,typename T3>
__global__ void copy_directions(T3* target_directions, const entryStructre_correlation<T,T3> *globalMem_correlation, 
        ub32 elements_num, bool is3D, bool is_mean)
{
    ub32 th_idx = threadIdx.x + blockIdx.x * blockDim.x;
    const entryStructre_correlation<T,T3> *current_element;

    // better in double to prevent data loss in middle calculations
    double3 dir_l;
    double mult_factor;

    if(th_idx < elements_num)
    {
        current_element = globalMem_correlation + th_idx;
        if(is_mean)
        {
            if(is3D)
            {
                mult_factor = (current_element->pixel_1.illumination_direction.z > 0.0 ? 1.0 : -1.0);
                dir_l.x = (current_element->pixel_1.illumination_direction.x + 
                        current_element->pixel_2.illumination_direction.x) / 2.0;
                dir_l.y = (current_element->pixel_1.illumination_direction.y + 
                        current_element->pixel_2.illumination_direction.y) / 2.0;
                dir_l.z = mult_factor * sqrt(1.0 - dir_l.x * dir_l.x - dir_l.y * dir_l.y);

                target_directions[th_idx].x = (T) dir_l.x;
                target_directions[th_idx].y = (T) dir_l.y;
                target_directions[th_idx].z = (T) dir_l.z;
            }
            else
            {
                mult_factor = (current_element->pixel_1.illumination_direction.y > 0.0 ? 1.0 : -1.0);
                dir_l.x = (current_element->pixel_1.illumination_direction.x + 
                        current_element->pixel_2.illumination_direction.x) / 2.0;
                dir_l.y = mult_factor * sqrt(1.0 - dir_l.x * dir_l.x);

                target_directions[th_idx].x = (T) dir_l.x;
                target_directions[th_idx].y = (T) dir_l.y;
            }
            
        }
        else
        {
            target_directions[th_idx] = current_element->pixel_1.illumination_direction;
            target_directions[th_idx + elements_num] = current_element->pixel_2.illumination_direction;
        }
    }
}

template<typename T,typename T3>
__global__ void copy_directions(T3* target_directions, const entryStructre_field<T,T3> *globalMem_field, ub32 elements_num)
{
    ub32 th_idx = threadIdx.x + blockIdx.x * blockDim.x;

    if(th_idx < elements_num)
    {
        target_directions[th_idx] = (globalMem_field + th_idx)->pixel.illumination_direction;
    }
}

template<typename T,typename T2,typename T3>
importanceSampling<T,T3>* sampling_preprocess(const input_IS<T,T2> *iS, const void *gpuGlobalData,
        ub32 is_correlation, bool is3D, ub32 elements_number)
{
    importanceSampling<T,T3> IS_struct, *gpu_IS_struct;

    bool is_mean = (bool) iS->is_mean;

    if(iS->direction_type == 1)
    {
        // Random sampling
        gpu_IS_struct = 0;
    }
    else if(iS->direction_type == 2)
    {
        // 2: Tabulated I.S. with f0
        IS_struct.f0_sampling = (tabulatedScatter<T>*) init_tabulated_scatterer<T,T2>(iS->tabulated_entries, iS->tabulated_values, is3D);
    }
    else if(iS->direction_type == 3)
    {
        // 3: HG I.S. with g0
        IS_struct.g0_sampling = (HGscatter<T>*) init_hg_scatterer<T>(iS->g0);
    }

    if(iS->direction_type != 1)
    {
        if(is_correlation)
        {
            const entryStructre_correlation<T,T3> *globalMem_correlation = (const entryStructre_correlation<T,T3> *) gpuGlobalData;
            IS_struct.directions_num = (is_mean ? 1 : 2) * elements_number;
            cudaMalloc(&IS_struct.directions_list, IS_struct.directions_num * sizeof(T3));
            copy_directions<T,T3><<<(elements_number - 1)/THREADS_NUM + 1,THREADS_NUM>>>(
                     IS_struct.directions_list, globalMem_correlation, elements_number, is3D, is_mean);
        }
        else
        {
            const entryStructre_field<T,T3> *globalMem_field = (const entryStructre_field<T,T3> *) gpuGlobalData;
            IS_struct.directions_num = elements_number;
            cudaMalloc(&IS_struct.directions_list, IS_struct.directions_num * sizeof(T3));
            copy_directions<T,T3><<<(elements_number - 1)/THREADS_NUM + 1,THREADS_NUM>>>(
                     IS_struct.directions_list, globalMem_field, elements_number);
        }

        cudaMalloc(&gpu_IS_struct, sizeof(importanceSampling<T,T3>));
        cudaMemcpy(gpu_IS_struct, &IS_struct, sizeof(importanceSampling<T,T3>), cudaMemcpyHostToDevice);
    }

    return gpu_IS_struct;
}

template<typename T, typename T3>
void sampling_preprocess_free(importanceSampling<T,T3> *IS_struct)
{
    if(IS_struct != 0)
    {
        importanceSampling<T,T3>* cpu_IS = (importanceSampling<T,T3>*) malloc(sizeof(importanceSampling<T,T3>));
        cudaMemcpy(cpu_IS,IS_struct, sizeof(importanceSampling<T,T3>),cudaMemcpyDeviceToHost);

        cudaFree(cpu_IS->g0_sampling);
        cudaFree(cpu_IS->directions_list);

        cudaFree(IS_struct);
        free(cpu_IS);
    }
}

template importanceSampling<double,double3>* sampling_preprocess<double,double2,double3>(const input_IS<double,double2> *iS, const void *gpuGlobalData,
        ub32 is_correlation, bool is3D, ub32 elements_number);
template importanceSampling<float,float3>* sampling_preprocess<float,float2,float3>(const input_IS<float,float2> *iS, const void *gpuGlobalData,
        ub32 is_correlation, bool is3D, ub32 elements_number);

template void sampling_preprocess_free<double,double3>(importanceSampling<double,double3> *IS_struct);
template void sampling_preprocess_free<float,float3>(importanceSampling<float,float3> *IS_struct);
