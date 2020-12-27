#include "refocus_header.h"
#include "refocus_scattering.h"
#include "cuMath.cu"

//debug
// #include<stdio.h>

__constant__ constStructre<double> constStr_scat_single;

template<typename T, typename T2, typename T3, typename scatteringType>
__global__ void singleScattering_kernel(T2* single_scattering_matrix, const T3* x0, 
	const scatteringType* scattering_struct, T k, const ff_direction<T3> *l_dir, const ff_direction<T3> *v_dir, ub32 directions_num)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_single;
    ub32 entry_num = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(entry_num < (directions_num * directions_num))
    {
        ub32 l_idx = entry_num % directions_num;
        ub32 v_idx = entry_num / directions_num;

        const ff_direction<T3> *il_dir = l_dir + l_idx;
        const ff_direction<T3> *iv_dir = v_dir + v_idx;

        if(il_dir->is_active && iv_dir->is_active)
        {
            T2 e;
            T bd1, bd2, bd3, dz;

            bd1 = abs(il_dir->dir_attenuation.x) < 1e-8 ? INFINITY : il_dir->dir_attenuation.x >= 0 ? (x0->x - curr_constStr->box_min[0]) / ( il_dir->dir_attenuation.x) : (x0->x - curr_constStr->box_max[0]) / (il_dir->dir_attenuation.x);
            bd2 = abs(il_dir->dir_attenuation.y) < 1e-8 ? INFINITY : il_dir->dir_attenuation.y >= 0 ? (x0->y - curr_constStr->box_min[1]) / ( il_dir->dir_attenuation.y) : (x0->y - curr_constStr->box_max[1]) / (il_dir->dir_attenuation.y); 
            bd3 = abs(il_dir->dir_attenuation.z) < 1e-8 ? INFINITY : il_dir->dir_attenuation.z >= 0 ? (x0->z - curr_constStr->box_min[2]) / ( il_dir->dir_attenuation.z) : (x0->z - curr_constStr->box_max[2]) / (il_dir->dir_attenuation.z);
            dz = fmin(fmin(bd1, bd2), bd3);

            bd1 = abs(iv_dir->dir_attenuation.x) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.x < 0 ? (curr_constStr->box_min[0] - x0->x) / ( iv_dir->dir_attenuation.x) : (curr_constStr->box_max[0] - x0->x) / (iv_dir->dir_attenuation.x); 
            bd2 = abs(iv_dir->dir_attenuation.y) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.y < 0 ? (curr_constStr->box_min[1] - x0->y) / ( iv_dir->dir_attenuation.y) : (curr_constStr->box_max[1] - x0->y) / (iv_dir->dir_attenuation.y); 
            bd3 = abs(iv_dir->dir_attenuation.z) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.z < 0 ? (curr_constStr->box_min[2] - x0->z) / ( iv_dir->dir_attenuation.z) : (curr_constStr->box_max[2] - x0->z) / (iv_dir->dir_attenuation.z);
            dz += fmin(fmin(bd1, bd2), bd3);

            e.x = -curr_constStr->sigt * dz;

            e.y = k * 
              (x0->x * (il_dir->dir.x - iv_dir->dir.x) +
               x0->y * (il_dir->dir.y - iv_dir->dir.y) +
               x0->z * (il_dir->dir.z - iv_dir->dir.z) );


            e = complexExponent(e);

            T s_lv = scattering_contribution<T,T3,3,true>(scattering_struct, &il_dir->dir, &iv_dir->dir);

            *(single_scattering_matrix + entry_num) = s_lv * e;
        }
        else
        {
            *(single_scattering_matrix + entry_num) = buildComplex((T) 0.0, (T) 0.0);
        }
    }
}

template<typename T, typename T2, typename T3>
void singleScattering(T2 *us, const T3 *x0, ff_refocus_struct<T2,T3>* ff_refocus, const T2* constPath, ub32 is_correlation, ub32 total_elements,
	T k, ub32 scattering_paramenter, const void* scatteringStruct)
{
	ub32 total_directions = ff_refocus->num_directions;
    // DEBUG

//     FILE *fptr;
// 
//     fptr = fopen("C:\\Users\\chen.bar\\Documents\\GitHub\\gaussianBeam-field\\mex\\debug_2.txt","w");
// {   
//     double2 A[12] =
//     {
//         {0.1,0.1}, {0.2,0.2}, {0.3,0.3}, {0.4,0.4},
//         {0.5,0.5}, {0.6,0.6}, {0.7,0.7}, {0.8,0.8},
//         {0.9,0.9}, {1.0,1.0}, {1.1,1.1}, {1.2,1.2}
//     };
// 
//     double2 B[20] =
//     {
//         {10.0,5.00}, {10.1,5.05}, {10.2,5.10}, {10.3,5.15}, {10.4,5.20},
//         {10.5,5.25}, {10.6,5.30}, {10.7,5.35}, {10.8,5.40}, {10.9,5.45},
//         {11.0,5.50}, {11.1,5.55}, {11.2,5.60}, {11.3,5.65}, {11.4,5.70},
//         {11.5,5.75}, {11.6,5.80}, {11.7,5.85}, {11.8,5.90}, {11.9,5.95}
//     };
// 
//     double2 C[15];
// 
//     double2 *gpu_a, *gpu_b, *gpu_c;
// 
//     cudaMalloc(&gpu_a, sizeof(double2) * 12);
//     cudaMalloc(&gpu_b, sizeof(double2) * 20);
//     cudaMalloc(&gpu_c, sizeof(double2) * 15);
// 
//     cudaMemcpy(gpu_a , A , sizeof(double2) * 12, cudaMemcpyHostToDevice);
//     cudaMemcpy(gpu_b , B , sizeof(double2) * 20, cudaMemcpyHostToDevice);
// 
//     matrixMult(ff_refocus->cublas_handle, 3, 4, 5, gpu_a, gpu_b, gpu_c);
//     cudaMemcpy(C , gpu_c , sizeof(double2) * 15, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"C_double = [");
//     for(ub32 iter1 = 0; iter1 < 3 ; iter1++)
//     {
//         for(ub32 iter2 = 0; iter2 < 5 ; iter2++)
//         {
//             fprintf(fptr,"%f + %fi ",C[iter2 + 5 * iter1].x,C[iter2 + 5 * iter1].y);
//         }
//         fprintf(fptr,"; \n");
//     }
//     fprintf(fptr,"];\n\n");
//     cudaFree(gpu_a);
//     cudaFree(gpu_b);
//     cudaFree(gpu_c);
// }
// 
// {   
//     float2 A[12] =
//     {
//         {0.1f,0.1f}, {0.2f,0.2f}, {0.3f,0.3f}, {0.4f,0.4f},
//         {0.5f,0.5f}, {0.6f,0.6f}, {0.7f,0.7f}, {0.8f,0.8f},
//         {0.9f,0.9f}, {1.0f,1.0f}, {1.1f,1.1f}, {1.2f,1.2f}
//     };
// 
//     float2 B[20] =
//     {
//         {10.0f,5.00f}, {10.1f,5.05f}, {10.2f,5.10f}, {10.3f,5.15f}, {10.4f,5.20f},
//         {10.5f,5.25f}, {10.6f,5.30f}, {10.7f,5.35f}, {10.8f,5.40f}, {10.9f,5.45f},
//         {11.0f,5.50f}, {11.1f,5.55f}, {11.2f,5.60f}, {11.3f,5.65f}, {11.4f,5.70f},
//         {11.5f,5.75f}, {11.6f,5.80f}, {11.7f,5.85f}, {11.8f,5.90f}, {11.9f,5.95f}
//     };
// 
//     float2 C[15];
// 
//     float2 *gpu_a, *gpu_b, *gpu_c;
// 
//     cudaMalloc(&gpu_a, sizeof(float2) * 12);
//     cudaMalloc(&gpu_b, sizeof(float2) * 20);
//     cudaMalloc(&gpu_c, sizeof(float2) * 15);
// 
//     cudaMemcpy(gpu_a , A , sizeof(float2) * 12, cudaMemcpyHostToDevice);
//     cudaMemcpy(gpu_b , B , sizeof(float2) * 20, cudaMemcpyHostToDevice);
// 
//     matrixMult(ff_refocus->cublas_handle, 3, 4, 5, gpu_a, gpu_b, gpu_c);
//     cudaMemcpy(C , gpu_c , sizeof(float2) * 15, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"C_single = [");
//     for(ub32 iter1 = 0; iter1 < 3 ; iter1++)
//     {
//         for(ub32 iter2 = 0; iter2 < 5 ; iter2++)
//         {
//             fprintf(fptr,"%f + %fi ",C[iter2 + 5 * iter1].x,C[iter2 + 5 * iter1].y);
//         }
//         fprintf(fptr,"; \n");
//     }
//     fprintf(fptr,"];\n\n");
//     cudaFree(gpu_a);
//     cudaFree(gpu_b);
//     cudaFree(gpu_c);
// }
// 
//     fclose(fptr);

    // END DEBUG
            
// DEBUG

//     FILE *fptr;
// 
//     fptr = fopen("C:\\Users\\chen.bar\\Documents\\GitHub\\gaussianBeam-field\\mex\\debug_2.txt","w");
// {   
//     double2 A[6] =
//     {
//         {0.1,0.1}, {0.2,0.2}, {0.3,0.3}, {0.4,0.4},
//         {0.5,0.5}, {0.6,0.6}
//     };
// 
//     double2 B[6] =
//     {
//         {10.0,5.00}, {10.1,5.05}, {10.2,5.10}, {10.3,5.15}, {10.4,5.20},
//         {10.5,5.25}
//     };
// 
//     double2 *gpu_a, *gpu_b;
// 
//     cudaMalloc(&gpu_a, sizeof(double2) * 6);
//     cudaMalloc(&gpu_b, sizeof(double2) * 6);
// 
//     cudaMemcpy(gpu_a , A , sizeof(double2) * 6, cudaMemcpyHostToDevice);
//     cudaMemcpy(gpu_b , B , sizeof(double2) * 6, cudaMemcpyHostToDevice);
// 
//     matrixAdd(ff_refocus->cublas_handle, 6, gpu_a, gpu_b);
//     cudaMemcpy(A , gpu_a , sizeof(double2) * 6, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"A_double = [");
//     for(ub32 iter1 = 0; iter1 < 6 ; iter1++)
//     {
//         fprintf(fptr,"%f + %fi ",A[iter1].x,A[iter1].y);
//     }
//     fprintf(fptr,"];\n\n");
//     cudaFree(gpu_a);
//     cudaFree(gpu_b);
// }
// 
// {   
//     float2 A[5] =
//     {
//         {0.1f,0.1f}, {0.2f,0.2f}, {0.3f,0.3f}, {0.4f,0.4f},
//         {0.5f,0.5f}
//     };
// 
//     float2 B[5] =
//     {
//         {10.0f,5.00f}, {10.1f,5.05f}, {10.2f,5.10f}, {10.3f,5.15f}, {10.4f,5.20f}
//     };
// 
//     float2 *gpu_a, *gpu_b;
// 
//     cudaMalloc(&gpu_a, sizeof(float2) * 5);
//     cudaMalloc(&gpu_b, sizeof(float2) * 5);
// 
//     cudaMemcpy(gpu_a , A , sizeof(float2) * 5, cudaMemcpyHostToDevice);
//     cudaMemcpy(gpu_b , B , sizeof(float2) * 5, cudaMemcpyHostToDevice);
// 
//     matrixAdd(ff_refocus->cublas_handle, 5, gpu_a, gpu_b);
//     cudaMemcpy(A , gpu_a , sizeof(float2) * 5, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"A_single = [");
//     for(ub32 iter1 = 0; iter1 < 5 ; iter1++)
//     {
//         fprintf(fptr,"%f + %fi ",A[iter1].x,A[iter1].y);
//     }
//     fprintf(fptr,"];\n\n");
//     cudaFree(gpu_a);
//     cudaFree(gpu_b);
// }
// 
//     fclose(fptr);

    // END DEBUG
            
// DEBUG
    
//     FILE *fptr;
// 
//     fptr = fopen("C:\\Users\\chen.bar\\Documents\\GitHub\\gaussianBeam-field\\mex\\debug_2.txt","w");
//     T2* wv_cpu, *wl_cpu;
// 
//     wv_cpu = (T2*) malloc(sizeof(T2) * total_elements * total_directions);
//     wl_cpu = (T2*) malloc(sizeof(T2) * total_elements * total_directions);
// 
//     cudaMemcpy(wv_cpu , ff_refocus->w_v , sizeof(T2) * total_elements * total_directions, cudaMemcpyDeviceToHost);
//     cudaMemcpy(wl_cpu , ff_refocus->w_l , sizeof(T2) * total_elements * total_directions, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"wv = [");
//     for(ub32 iter1 = 0; iter1 < total_elements ; iter1++)
//     {
//         for(ub32 iter2 = 0; iter2 < total_directions ; iter2++)
//         {
//             fprintf(fptr,"%e + %ei ",wv_cpu[iter2 + total_directions * iter1].x,wv_cpu[iter2 + total_directions * iter1].y);
//         }
//         fprintf(fptr,"; \n");
//     }
//     fprintf(fptr,"];\n\n");
// 
//     fprintf(fptr,"wl = [");
//     for(ub32 iter1 = 0; iter1 < total_elements ; iter1++)
//     {
//         for(ub32 iter2 = 0; iter2 < total_directions ; iter2++)
//         {
//             fprintf(fptr,"%e + %ei ",wl_cpu[iter2 + total_directions * iter1].x,wl_cpu[iter2 + total_directions * iter1].y);
//         }
//         fprintf(fptr,"; \n");
//     }
//     fprintf(fptr,"];\n\n");
// 
// 
//     free(wv_cpu);
//     free(wl_cpu);
// 
//     fclose(fptr);

//     END DEBUG
            
    ub32 blocks_num = (total_directions * total_directions - 1) / THREADS_NUM + 1;
    ub32 elementwise_blocks_num = (total_elements * total_directions - 1) / THREADS_NUM + 1;

    // build the ff matrix
    if(scattering_paramenter == 1)
    {
        singleScattering_kernel<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>
            (ff_refocus->single_scattering_ff,x0,(const randomScatter*) scatteringStruct,k,
             ff_refocus->d_l,ff_refocus->d_v,total_directions);
    }

    if(scattering_paramenter == 2)
    {
        singleScattering_kernel<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>
            (ff_refocus->single_scattering_ff,x0,(tabulatedScatter<T>*) scatteringStruct,k,
             ff_refocus->d_l,ff_refocus->d_v,total_directions);
    }

    if(scattering_paramenter == 3)
    {
        singleScattering_kernel<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>
            (ff_refocus->single_scattering_ff,x0,(HGscatter<T>*) scatteringStruct,k,
             ff_refocus->d_l,ff_refocus->d_v,total_directions);
    }

// DEBUG
    
//     FILE *fptr;
// 
//     fptr = fopen("C:\\Users\\chen.bar\\Documents\\GitHub\\gaussianBeam-field\\mex\\debug_2.txt","w");
//     T2* ff_mat_cpu;
// 
//     ff_mat_cpu = (T2*) malloc(sizeof(T2) * total_directions * total_directions);
// 
//     cudaMemcpy(ff_mat_cpu , ff_refocus->single_scattering_ff , sizeof(T2) * total_directions * total_directions, cudaMemcpyDeviceToHost);
// 
//     fprintf(fptr,"ff_mat = [");
//     for(ub32 iter1 = 0; iter1 < total_directions ; iter1++)
//     {
//         for(ub32 iter2 = 0; iter2 < total_directions ; iter2++)
//         {
//             fprintf(fptr,"%e + %ei ",ff_mat_cpu[iter2 + total_directions * iter1].x,ff_mat_cpu[iter2 + total_directions * iter1].y);
//         }
//         fprintf(fptr,"; \n");
//     }
//     fprintf(fptr,"];\n\n");
// 
//     free(ff_mat_cpu);
// 
//     fclose(fptr);

//     END DEBUG

    // wv * scattering matrix          
	matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, total_directions, 
        ff_refocus->w_v, ff_refocus->single_scattering_ff, ff_refocus->mult_w_v_single_scattering);

    // (wv * scattering matrix) .* wl
    elementwise_mult<T2><<<elementwise_blocks_num,THREADS_NUM>>>(ff_refocus->mult_w_v_single_scattering, ff_refocus->w_l, total_elements * total_directions);

    // sum((wv * scattering matrix) .* wl,2)
    matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
        ff_refocus->mult_w_v_single_scattering, ff_refocus->ones_vector, ff_refocus->mult_res);

    if(is_correlation)
    {
        // build the ff matrix
        if(scattering_paramenter == 1)
        {
            singleScattering_kernel<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>
                (ff_refocus->single_scattering_ff,x0,(const randomScatter*) scatteringStruct,k,
                 ff_refocus->d_l_2,ff_refocus->d_v_2,total_directions);
        }

        if(scattering_paramenter == 2)
        {
            singleScattering_kernel<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>
                (ff_refocus->single_scattering_ff,x0,(tabulatedScatter<T>*) scatteringStruct,k,
                 ff_refocus->d_l_2,ff_refocus->d_v_2,total_directions);
        }

        if(scattering_paramenter == 3)
        {
            singleScattering_kernel<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>
                (ff_refocus->single_scattering_ff,x0,(HGscatter<T>*) scatteringStruct,k,
                 ff_refocus->d_l_2,ff_refocus->d_v_2,total_directions);
        }

        // wv_2 * scattering matrix
        matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, total_directions, 
            ff_refocus->w_v_2, ff_refocus->single_scattering_ff, ff_refocus->mult_w_v_single_scattering);


        // (wv_2 * scattering matrix) .* wl_2
        elementwise_mult<T2><<<elementwise_blocks_num,THREADS_NUM>>>(ff_refocus->mult_w_v_single_scattering, ff_refocus->w_l_2, total_elements * total_directions);

        // sum((wv * scattering matrix) .* wl_2,2)
        matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
            ff_refocus->mult_w_v_single_scattering, ff_refocus->ones_vector, ff_refocus->mult_res_2);

        ub32 mult_blocks_num = (total_elements - 1) / THREADS_NUM + 1;
        elementwise_conj_mult<T2><<<mult_blocks_num,THREADS_NUM>>>(ff_refocus->mult_res, ff_refocus->mult_res_2, total_elements);
    }

    matrixAdd(ff_refocus->cublas_handle_device_alpha, total_elements,
        constPath, us, ff_refocus->mult_res);
}

template<typename T>
void singleScattering_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_scat_single, constMem, sizeof(constStructre<T>));
}

template void singleScattering_setConstMem<double>(constStructre<double> *constMem);
template void singleScattering_setConstMem<float>(constStructre<float> *constMem);

template void singleScattering<double,double2,double3>(double2 *us, const double3 *x0,
    ff_refocus_struct<double2,double3>* ff_refocus, const double2* constPath,
    ub32 is_correlation, ub32 total_elements,
	double k, ub32 scattering_paramenter, const void* scatteringStruct);

template void singleScattering<float,float2,float3>(float2 *us, const float3 *x0,
    ff_refocus_struct<float2,float3>* ff_refocus, const float2* constPath,
    ub32 is_correlation, ub32 total_elements,
	float k, ub32 scattering_paramenter, const void* scatteringStruct);
