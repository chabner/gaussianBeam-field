#include "ff_header.h"
#include "ff_scattering.h"
#include "cuMath.cu"

__constant__ constStructre<double> constStr_scat_mult;
__constant__ ub32 elemNum_scat_mult;
__constant__ ub32 is_inside_elemes;

template<typename T, typename T2, typename T3, typename scatteringType, ub32 dims, bool is_cbs>
__device__ T2 el_times_ev(const T3 x0, const T3 w0, T3 xb, T3 wb, T path_length,
        const pixel_entry<T,T3>* pixelData, const scatteringType* scatteringStruct)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_mult;

    T2 e;

    T bd1, bd2, bd3, dz;
    T sct_l, sct_v;

    T2 res;

    bd1 = abs(pixelData->illumination_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.x >= 0 ? (x0.x - curr_constStr->box_min[0]) / ( pixelData->illumination_attenuation_direction.x) : (x0.x - curr_constStr->box_max[0]) / (pixelData->illumination_attenuation_direction.x);
    bd2 = abs(pixelData->illumination_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.y >= 0 ? (x0.y - curr_constStr->box_min[1]) / ( pixelData->illumination_attenuation_direction.y) : (x0.y - curr_constStr->box_max[1]) / (pixelData->illumination_attenuation_direction.y); 
    if(dims == 3)
    {
        bd3 = abs(pixelData->illumination_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.z >= 0 ? (x0.z - curr_constStr->box_min[2]) / ( pixelData->illumination_attenuation_direction.z) : (x0.z - curr_constStr->box_max[2]) / (pixelData->illumination_attenuation_direction.z);
        dz = fmin(fmin(bd1, bd2), bd3);
    }
    if(dims == 2)
    {
        dz = fmin(bd1, bd2);
    }

    bd1 = abs(pixelData->view_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.x < 0 ? (curr_constStr->box_min[0] - xb.x) / ( pixelData->view_attenuation_direction.x) : (curr_constStr->box_max[0] - xb.x) / (pixelData->view_attenuation_direction.x); 
    bd2 = abs(pixelData->view_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.y < 0 ? (curr_constStr->box_min[1] - xb.y) / ( pixelData->view_attenuation_direction.y) : (curr_constStr->box_max[1] - xb.y) / (pixelData->view_attenuation_direction.y); 
    if(dims == 3)
    {
        T bd3 = abs(pixelData->view_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.z < 0 ? (curr_constStr->box_min[2] - xb.z) / ( pixelData->view_attenuation_direction.z) : (curr_constStr->box_max[2] - xb.z) / (pixelData->view_attenuation_direction.z);
        dz += fmin(fmin(bd1, bd2), bd3);
    }
    if(dims == 2)
    {
        dz += fmin(bd1, bd2);
    }

    e.x = -curr_constStr->sigt * dz;
    if(dims == 3)
    {
        e.y = pixelData->k * 
                  (x0.x * pixelData->illumination_direction.x - xb.x * pixelData->view_direction.x +
                   x0.y * pixelData->illumination_direction.y - xb.y * pixelData->view_direction.y +
                   x0.z * pixelData->illumination_direction.z - xb.z * pixelData->view_direction.z + 
                   path_length);
    }
    if(dims == 2)
    {
        e.y = pixelData->k * 
                  (x0.x * pixelData->illumination_direction.x - xb.x * pixelData->view_direction.x +
                   x0.y * pixelData->illumination_direction.y - xb.y * pixelData->view_direction.y +
                   path_length);
    }

    sct_l = scattering_contribution<T,T3,dims,true>(scatteringStruct,&pixelData->illumination_direction,&w0);
    sct_v = scattering_contribution<T,T3,dims,true>(scatteringStruct,&pixelData->view_direction,&wb);

    res = (sct_l * sct_v) * complexExponent(e);

    if(is_cbs)
    {
        bd1 = abs(pixelData->illumination_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.x >= 0 ? (xb.x - curr_constStr->box_min[0]) / ( pixelData->illumination_attenuation_direction.x) : (xb.x - curr_constStr->box_max[0]) / (pixelData->illumination_attenuation_direction.x);
        bd2 = abs(pixelData->illumination_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.y >= 0 ? (xb.y - curr_constStr->box_min[1]) / ( pixelData->illumination_attenuation_direction.y) : (xb.y - curr_constStr->box_max[1]) / (pixelData->illumination_attenuation_direction.y); 
        if(dims == 3)
        {
            bd3 = abs(pixelData->illumination_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.z >= 0 ? (xb.z - curr_constStr->box_min[2]) / ( pixelData->illumination_attenuation_direction.z) : (xb.z - curr_constStr->box_max[2]) / (pixelData->illumination_attenuation_direction.z);
            dz = fmin(fmin(bd1, bd2), bd3);
        }
        if(dims == 2)
        {
            dz = fmin(bd1, bd2);
        }

        bd1 = abs(pixelData->view_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.x < 0 ? (curr_constStr->box_min[0] - x0.x) / ( pixelData->view_attenuation_direction.x) : (curr_constStr->box_max[0] - x0.x) / (pixelData->view_attenuation_direction.x); 
        bd2 = abs(pixelData->view_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.y < 0 ? (curr_constStr->box_min[1] - x0.y) / ( pixelData->view_attenuation_direction.y) : (curr_constStr->box_max[1] - x0.y) / (pixelData->view_attenuation_direction.y); 
        if(dims == 3)
        {
            T bd3 = abs(pixelData->view_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.z < 0 ? (curr_constStr->box_min[2] - x0.z) / ( pixelData->view_attenuation_direction.z) : (curr_constStr->box_max[2] - x0.z) / (pixelData->view_attenuation_direction.z);
            dz += fmin(fmin(bd1, bd2), bd3);
        }
        if(dims == 2)
        {
            dz += fmin(bd1, bd2);
        }

        e.x = -curr_constStr->sigt * dz;
        if(dims == 3)
        {
            e.y = pixelData->k * 
                      (xb.x * pixelData->illumination_direction.x - x0.x * pixelData->view_direction.x +
                       xb.y * pixelData->illumination_direction.y - x0.y * pixelData->view_direction.y +
                       xb.z * pixelData->illumination_direction.z - x0.z * pixelData->view_direction.z + 
                       path_length);
        }
        if(dims == 2)
        {
            e.y = pixelData->k * 
                      (xb.x * pixelData->illumination_direction.x - x0.x * pixelData->view_direction.x +
                       xb.y * pixelData->illumination_direction.y - x0.y * pixelData->view_direction.y +
                       path_length);
        }

        T3 wb2, w02;
        wb2.x = ((T)-1.0) * wb.x;
        wb2.y = ((T)-1.0) * wb.y;
        wb2.z = ((T)-1.0) * wb.z;

        w02.x = ((T)-1.0) * w0.x;
        w02.y = ((T)-1.0) * w0.y;
        w02.z = ((T)-1.0) * w0.z;

        sct_l = scattering_contribution<T,T3,dims,true>(scatteringStruct,&pixelData->illumination_direction,&wb2);
        sct_v = scattering_contribution<T,T3,dims,true>(scatteringStruct,&pixelData->view_direction,&w02);

        res = (res + (sct_l * sct_v) * complexExponent(e)) / ((T) (sqrt(2.0)));
    }

    return res;
}

template<typename T, typename T2, typename T3, typename scatteringType, ub32 dims, bool is_cbs>
__device__  void scatteringLoop(ub32 current_path, ub32 current_pixel_offset,
        const T3* x0, const T3* xb, const T3* w0, const T3* wb,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_correlation<T,T3> *pixelData,
        const T2* constPath, const T* path_length, const scatteringType* scatteringStruct)
{      
    T2 res   = el_times_ev<T,T2,T3,scatteringType,dims,is_cbs>(
            x0[current_path],
            w0[current_path],
            xb[current_path],
            wb[current_path],
            path_length[current_path],
            &(pixelData + current_pixel_offset)->pixel_1,
            scatteringStruct);

    T2 res_2 = el_times_ev<T,T2,T3,scatteringType,dims,is_cbs>(
            x0[current_path],
            w0[current_path],
            xb[current_path],
            wb[current_path],
            path_length[current_path],
            &(pixelData + current_pixel_offset)->pixel_2,
            scatteringStruct);
    T2 u_mult = conjMult(res,res_2);

    // Each thread copies the sum of all mixture to the shared memory
    u_mult = u_mult * constPath[current_path];
    u_res_x[threadIdx.x] = u_mult.x;
    u_res_y[threadIdx.x] = u_mult.y;
}

template<typename T, typename T2, typename T3, typename scatteringType, ub32 dims, bool is_cbs>
__device__ void scatteringLoop(ub32 current_path, ub32 current_pixel_offset,
        const T3* x0, const T3* xb, const T3* w0, const T3* wb,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_field<T,T3> *pixelData,
        const T2* constPath, const T* path_length, const scatteringType* scatteringStruct)
{
    T2 res = el_times_ev<T,T2,T3,scatteringType,dims,is_cbs>(
        x0[current_path],
        w0[current_path],
        xb[current_path],
        wb[current_path],
        path_length[current_path],
        &(pixelData + current_pixel_offset)->pixel,
        scatteringStruct);
    

    res = res * constPath[current_path];
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;
}

template<typename T, typename T2, typename T3, typename correlationType, typename scatteringType, ub32 dims, bool is_cbs, ub32 blocksNum>
__launch_bounds__(THREADS_NUM)
__global__ void multipleScattering_kernel(T2* u,
        const T3* x0, const T3* xb, const T3* w0, const T3* wb,
        const correlationType* dataIn, const T2* constPath, const scatteringType* scatteringStruct,
        const ub32* is_inside, const T* path_length)
{
    __shared__ volatile T u_res_x[THREADS_NUM];
    __shared__ volatile T u_res_y[THREADS_NUM];
    __shared__ correlationType pixelData[blocksNum];

    ub32 firstPixel = blocksNum * blockIdx.x;
    ub32 lastPixel = firstPixel + (blocksNum - 1);
    if(lastPixel >= elemNum_scat_mult) {lastPixel = elemNum_scat_mult - 1;}
    ub32 pixelsPerBlock = lastPixel - firstPixel + 1;

    // copy the relevant data to shared memory
    int* iDest = (int*)pixelData;
    const int* iSrc = (const int*)(dataIn + firstPixel);

    for(unsigned int i = threadIdx.x; i<pixelsPerBlock*sizeof(correlationType)/sizeof(int); i+=blockDim.x)
        iDest[i] = iSrc[i];

    __syncthreads();

    // check which pixel and path the current thread is on
    ub32 current_pixel_offset = threadIdx.x % blocksNum;
    ub32 current_path_num = threadIdx.x / blocksNum;
    bool is_legal_pixel = (current_pixel_offset < pixelsPerBlock) & (current_path_num < is_inside_elemes);

    if(is_legal_pixel)
    {
        ub32 current_path = is_inside[current_path_num];
        scatteringLoop<T,T2,T3,scatteringType,dims,is_cbs>(current_path, current_pixel_offset, x0, xb, w0, wb, u_res_x, u_res_y, pixelData, constPath, path_length, scatteringStruct);
    }
    else
    {
        u_res_x[threadIdx.x] = (T) 0.0;
        u_res_y[threadIdx.x] = (T) 0.0;
    }

    __syncthreads();

    // Reduce sum of all shared memory to a single result
    if(blocksNum <= 512){if (threadIdx.x < 512){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 512];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 512];
    }  __syncthreads();}
    if(blocksNum <= 256){if (threadIdx.x < 256){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 256];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 256];
    }  __syncthreads();}
    if(blocksNum <= 128){if (threadIdx.x < 128){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 128];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 128];
    }  __syncthreads();}
    if(blocksNum <= 64) {if (threadIdx.x < 64 ){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 64 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 64 ];
    }  __syncthreads();}

    if(blocksNum <= 32)
    {
        if (threadIdx.x < 32 ) // warpReduce
        {
            u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 32];
            u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 32];

            if(blocksNum <= 16)
            {
                u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 16];
                u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 16];
            }

            if(blocksNum <= 8)
            {
                u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 8 ];
                u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 8 ];
            }

            if(blocksNum <= 4)
            {
                u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 4 ];
                u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 4 ];
            }

            if(blocksNum <= 2)
            {
                u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 2 ];
                u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 2 ];
            }

            if(blocksNum <= 1)
            {
                u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 1 ];
                u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 1 ];
            }
        }
    }

//     copy the final result to the global mem
    if (threadIdx.x < pixelsPerBlock) { u[firstPixel + threadIdx.x].x += u_res_x[threadIdx.x]; u[firstPixel + threadIdx.x].y += u_res_y[threadIdx.x]; }

}

template<typename T, typename T2, typename T3, typename correlationType, typename scatteringType, ub32 dims, bool is_cbs>
void ms_inner(T2 *u, const T3 *x0, const T3 *xb, const T3 *w0, const T3 *wb, const T2 *constPath, const pathsList<T>* pl, 
        ub32 total_elements, ub32 total_paths_number, const correlationType* globalMem, const scatteringType* scatteringStruct)
{
    if(total_paths_number > 512)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,1><<<total_elements, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 256)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,2><<<(total_elements - 1)/2 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 128)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,4><<<(total_elements - 1)/4 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 64)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,8><<<(total_elements - 1)/8 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 32)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,16><<<(total_elements - 1)/16 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 16)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,32><<<(total_elements - 1)/32 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else if(total_paths_number > 8)
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,64><<<(total_elements - 1)/64 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
    else
    {
        multipleScattering_kernel<T,T2,T3,correlationType,scatteringType,dims,is_cbs,128><<<(total_elements - 1)/128 + 1, THREADS_NUM>>>
                (u,x0,xb,w0,wb,globalMem,constPath,scatteringStruct,pl->paths_numbers_list,pl->path_length);
    }
}

template<typename T, typename T2, typename T3>
void multipleScattering(T2 *u, const T3 *x0, const T3 *xb, const T3 *w0, const T3 *wb, const T2 *constPath, const pathsList<T>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, ub32 dims, ub32 is_cbs, const void* globalMem, 
        ub32 scattering_paramenter, const void* scatteringStruct)
{
    cudaMemcpyToSymbol(is_inside_elemes, &total_paths_number, sizeof(ub32));

	if(is_correlation)
    { 
        if(scattering_paramenter == 1)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
            }
        }

        if(scattering_paramenter == 2)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
            }
        }

        if(scattering_paramenter == 3)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_correlation<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
            }
        }
    }
    else
    { 
        if(scattering_paramenter == 1)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,randomScatter,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,randomScatter,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,randomScatter,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,randomScatter,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const randomScatter*) scatteringStruct);
                }
            }
        }

        if(scattering_paramenter == 2)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const tabulatedScatter<T>*) scatteringStruct);
                }

            }
        }

        if(scattering_paramenter == 3)
        {
            if(dims == 2)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,2,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,2,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
            }
            if(dims == 3)
            {
                if(is_cbs)
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,3,true>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
                else
                {
                    ms_inner<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,3,false>(u,x0,xb,w0,wb,constPath,pl,
                        total_elements,total_paths_number,
                        (const entryStructre_field<T,T3>*) globalMem,(const HGscatter<T>*) scatteringStruct);
                }
            }
        }
    }
}

template<typename T>
void multipleScattering_setConstMem(constStructre<T> *constMem, ub32 number_of_elements)
{
    cudaMemcpyToSymbol(constStr_scat_mult, constMem, sizeof(constStructre<T>));
    cudaMemcpyToSymbol(elemNum_scat_mult, &number_of_elements, sizeof(ub32));
}

template void multipleScattering_setConstMem<double>(constStructre<double> *constMem, ub32 number_of_elements);
template void multipleScattering_setConstMem<float>(constStructre<float> *constMem, ub32 number_of_elements);

template void multipleScattering<double,double2,double3>(double2 *u, const double3 *x0, const double3 *xb, const double3 *w0,
        const double3 *wb, const double2 *constPath, const pathsList<double>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, ub32 dims, ub32 is_cbs, const void* globalMem,
        ub32 scattering_paramenter, const void* scatteringStruct);

template void multipleScattering<float,float2,float3>(float2 *u, const float3 *x0, const float3 *xb, const float3 *w0,
        const float3 *wb, const float2 *constPath, const pathsList<float>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, ub32 dims, ub32 is_cbs, const void* globalMem, 
        ub32 scattering_paramenter, const void* scatteringStruct);
        