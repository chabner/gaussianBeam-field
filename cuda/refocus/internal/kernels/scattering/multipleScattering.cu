#include "refocus_header.h"
#include "refocus_scattering.h"
#include "cuMath.cu"

__constant__ constStructre<double> constStr_scat_mult;

template<typename T, typename T2, typename T3, typename scatteringType>
__global__ void el_matrix(T2* el_matrix_res, const T3* x0, const T3* w0,
	const scatteringType* scatteringStruct, T k, const ff_direction<T3> *l_dir, ub32 directions_num)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_mult;
    ub32 entry_num = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(entry_num < directions_num)
    {
        const ff_direction<T3> *il_dir = l_dir + entry_num;

        if(il_dir->is_active)
        {
            T2 e;
            T bd1, bd2, bd3, dz, sct_l;

            bd1 = abs(il_dir->dir_attenuation.x) < 1e-8 ? INFINITY : il_dir->dir_attenuation.x >= 0 ? (x0->x - curr_constStr->box_min[0]) / ( il_dir->dir_attenuation.x) : (x0->x - curr_constStr->box_max[0]) / (il_dir->dir_attenuation.x);
            bd2 = abs(il_dir->dir_attenuation.y) < 1e-8 ? INFINITY : il_dir->dir_attenuation.y >= 0 ? (x0->y - curr_constStr->box_min[1]) / ( il_dir->dir_attenuation.y) : (x0->y - curr_constStr->box_max[1]) / (il_dir->dir_attenuation.y); 
            bd3 = abs(il_dir->dir_attenuation.z) < 1e-8 ? INFINITY : il_dir->dir_attenuation.z >= 0 ? (x0->z - curr_constStr->box_min[2]) / ( il_dir->dir_attenuation.z) : (x0->z - curr_constStr->box_max[2]) / (il_dir->dir_attenuation.z);
            dz = fmin(fmin(bd1, bd2), bd3);

            e.x = -curr_constStr->sigt * dz;
            e.y = k * 
              (x0->x * il_dir->dir.x +
               x0->y * il_dir->dir.y +
               x0->z * il_dir->dir.z);

            sct_l = scattering_contribution<T,T3,3,true>(scatteringStruct,&il_dir->dir,w0);

            el_matrix_res[entry_num] = sct_l * complexExponent(e);
        }
    }
}

template<typename T, typename T2, typename T3, typename scatteringType>
__global__ void ev_matrix(T2* ev_matrix_res, const T3* xb, const T3* wb,
	const scatteringType* scatteringStruct, T k, const ff_direction<T3> *v_dir, ub32 directions_num)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_mult;
    ub32 entry_num = blockIdx.x * blockDim.x + threadIdx.x;
    
    if(entry_num < directions_num)
    {
        const ff_direction<T3> *iv_dir = v_dir + entry_num;

        if(iv_dir->is_active)
        {
            T2 e;
            T bd1, bd2, bd3, dz, sct_v;

            bd1 = abs(iv_dir->dir_attenuation.x) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.x < 0 ? (curr_constStr->box_min[0] - xb->x) / ( iv_dir->dir_attenuation.x) : (curr_constStr->box_max[0] - xb->x) / (iv_dir->dir_attenuation.x); 
            bd2 = abs(iv_dir->dir_attenuation.y) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.y < 0 ? (curr_constStr->box_min[1] - xb->y) / ( iv_dir->dir_attenuation.y) : (curr_constStr->box_max[1] - xb->y) / (iv_dir->dir_attenuation.y); 
            bd3 = abs(iv_dir->dir_attenuation.z) < 1e-8 ? INFINITY : iv_dir->dir_attenuation.z < 0 ? (curr_constStr->box_min[2] - xb->z) / ( iv_dir->dir_attenuation.z) : (curr_constStr->box_max[2] - xb->z) / (iv_dir->dir_attenuation.z);
            dz = fmin(fmin(bd1, bd2), bd3);

            e.x = -curr_constStr->sigt * dz;
            e.y = -k * 
              (xb->x * iv_dir->dir.x +
               xb->y * iv_dir->dir.y +
               xb->z * iv_dir->dir.z);

            sct_v = scattering_contribution<T,T3,3,true>(scatteringStruct,&iv_dir->dir,wb);

            ev_matrix_res[entry_num] = sct_v * complexExponent(e);
        }
    }
}

template<typename T, typename T2, typename T3>
void multipleScattering(T2 *u, const ff_refocus_struct<T2,T3>* ff_refocus, T k,
                        ub32 scattering_paramenter, const void* scatteringStruct,
                        const T3 *x0, const T3* xb, const T3* w0, const T3* wb, const T2* constPath,
                        ub32 is_correlation, ub32 total_elements)
{
    ub32 total_directions = ff_refocus->num_directions;
    ub32 blocks_num = (total_directions - 1) / THREADS_NUM + 1;

    ub32 elementwise_blocks_num = (total_elements - 1) / THREADS_NUM + 1;

    // build el ff matrix
    if(scattering_paramenter == 1)
    {
        el_matrix<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_illumination_ff, x0, w0,
            (const randomScatter*) scatteringStruct, k, ff_refocus->d_l, total_directions);

        ev_matrix<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_view_ff, xb, wb,
            (const randomScatter*) scatteringStruct, k, ff_refocus->d_v, total_directions);
    }

    if(scattering_paramenter == 2)
    {
        el_matrix<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_illumination_ff, x0, w0,
            (const tabulatedScatter<T>*) scatteringStruct, k, ff_refocus->d_l, total_directions);

        ev_matrix<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_view_ff, xb, wb,
            (const tabulatedScatter<T>*) scatteringStruct, k, ff_refocus->d_v, total_directions);
    }

    if(scattering_paramenter == 3)
    {
        el_matrix<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_illumination_ff, x0, w0,
            (const HGscatter<T>*) scatteringStruct, k, ff_refocus->d_l, total_directions);

        ev_matrix<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>(
            ff_refocus->multiple_scattering_view_ff, xb, wb,
            (const HGscatter<T>*) scatteringStruct, k, ff_refocus->d_v, total_directions);
    }

    // wv * scattering matrix
	matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
        ff_refocus->w_v, ff_refocus->multiple_scattering_view_ff, ff_refocus->mult_w_v_multiple_scattering);

    // wl * scattering matrix
	matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
        ff_refocus->w_l, ff_refocus->multiple_scattering_illumination_ff, ff_refocus->mult_w_l_multiple_scattering);

    // (wv * scattering matrix) .* (wl * scattering matrix)
    elementwise_mult<T2><<<elementwise_blocks_num,THREADS_NUM>>>(ff_refocus->mult_w_v_multiple_scattering, ff_refocus->mult_w_l_multiple_scattering, total_elements);
    cudaMemcpy(ff_refocus->mult_res, ff_refocus->mult_w_v_multiple_scattering, total_elements * sizeof(T2), cudaMemcpyDeviceToDevice);

    if(is_correlation)
    {
        // build el ff matrix
        if(scattering_paramenter == 1)
        {
            el_matrix<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_illumination_ff, x0, w0,
                (const randomScatter*) scatteringStruct, k, ff_refocus->d_l_2, total_directions);

            ev_matrix<T,T2,T3,randomScatter><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_view_ff, xb, wb,
                (const randomScatter*) scatteringStruct, k, ff_refocus->d_v_2, total_directions);
        }

        if(scattering_paramenter == 2)
        {
            el_matrix<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_illumination_ff, x0, w0,
                (const tabulatedScatter<T>*) scatteringStruct, k, ff_refocus->d_l_2, total_directions);

            ev_matrix<T,T2,T3,tabulatedScatter<T>><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_view_ff, xb, wb,
                (const tabulatedScatter<T>*) scatteringStruct, k, ff_refocus->d_v_2, total_directions);
        }

        if(scattering_paramenter == 3)
        {
            el_matrix<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_illumination_ff, x0, w0,
                (const HGscatter<T>*) scatteringStruct, k, ff_refocus->d_l_2, total_directions);

            ev_matrix<T,T2,T3,HGscatter<T>><<<blocks_num, THREADS_NUM>>>(
                ff_refocus->multiple_scattering_view_ff, xb, wb,
                (const HGscatter<T>*) scatteringStruct, k, ff_refocus->d_v_2, total_directions);
        }

        // wv * scattering matrix
        matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
            ff_refocus->w_v_2, ff_refocus->multiple_scattering_view_ff, ff_refocus->mult_w_v_multiple_scattering);

        // wl * scattering matrix
        matrixMult(ff_refocus->cublas_handle, total_elements, total_directions, 1, 
            ff_refocus->w_l_2, ff_refocus->multiple_scattering_illumination_ff, ff_refocus->mult_w_l_multiple_scattering);

        // (wv * scattering matrix) .* (wl * scattering matrix)
        elementwise_mult<T2><<<elementwise_blocks_num,THREADS_NUM>>>(ff_refocus->mult_w_v_multiple_scattering, ff_refocus->mult_w_l_multiple_scattering, total_elements);
        cudaMemcpy(ff_refocus->mult_res_2, ff_refocus->mult_w_v_multiple_scattering, total_elements * sizeof(T2), cudaMemcpyDeviceToDevice);

        elementwise_conj_mult<T2><<<elementwise_blocks_num,THREADS_NUM>>>(ff_refocus->mult_res, ff_refocus->mult_res_2, total_elements);
    }

    matrixAdd(ff_refocus->cublas_handle_device_alpha, total_elements,
        constPath, u, ff_refocus->mult_res);
}

template<typename T>
void multipleScattering_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_scat_mult, constMem, sizeof(constStructre<T>));
}

template void multipleScattering<double,double2,double3>(double2 *u, const ff_refocus_struct<double2,double3>* ff_refocus, double k,
    ub32 scattering_paramenter, const void* scatteringStruct,
    const double3 *x0, const double3* xb, const double3* w0, const double3* wb, const double2* constPath,
    ub32 is_correlation, ub32 total_elements);

template void multipleScattering<float,float2,float3>(float2 *u, const ff_refocus_struct<float2,float3>* ff_refocus, float k,
    ub32 scattering_paramenter, const void* scatteringStruct,
    const float3 *x0, const float3* xb, const float3* w0, const float3* wb, const float2* constPath,
    ub32 is_correlation, ub32 total_elements);

template void multipleScattering_setConstMem<double>(constStructre<double> *constMem);
template void multipleScattering_setConstMem<float>(constStructre<float> *constMem);
