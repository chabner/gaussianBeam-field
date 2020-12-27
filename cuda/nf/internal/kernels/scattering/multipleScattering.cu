#include "nf_header.h"
#include "nf_scattering.h"
#include "cuMath.cu"

__constant__ constStructre<double> constStr_scat_mult;
__constant__ ub32 elemNum_scat_mult;
__constant__ ub32 is_inside_elemes;

template<typename T, typename T2, typename T3>
__device__ T2 el(const T3 x0, const T3 w0,
        const pixel_entry<T,T3>* __restrict__ pixelData)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_mult;
    T2 el = {0};
    T2 thL[3];
    T2 sqrtMu;

    T bd1, bd2, bd3, dz;
    unsigned char mixtureIdx;
    T gamma_s;
    T th_c, log_nu;

    // Distance to edge of sample
    bd1 = abs(pixelData->illuminationDir.x) < 1e-8 ? INFINITY : pixelData->illuminationDir.x >= 0 ? (x0.x - curr_constStr->box_min[0]) / ( pixelData->illuminationDir.x) : (x0.x - curr_constStr->box_max[0]) / (pixelData->illuminationDir.x);
    bd2 = abs(pixelData->illuminationDir.y) < 1e-8 ? INFINITY : pixelData->illuminationDir.y >= 0 ? (x0.y - curr_constStr->box_min[1]) / ( pixelData->illuminationDir.y) : (x0.y - curr_constStr->box_max[1]) / (pixelData->illuminationDir.y); 
    bd3 = abs(pixelData->illuminationDir.z) < 1e-8 ? INFINITY : pixelData->illuminationDir.z >= 0 ? (x0.z - curr_constStr->box_min[2]) / ( pixelData->illuminationDir.z) : (x0.z - curr_constStr->box_max[2]) / (pixelData->illuminationDir.z);

    dz  = fmin(fmin(bd1, bd2), bd3);

    // The complex throughput
    thL[0].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.x;
    thL[1].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.y;
    thL[2].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.z;

    thL[0].y = pixelData->k * (x0.x - pixelData->illuminationP.x);
    thL[1].y = pixelData->k * (x0.y - pixelData->illuminationP.y);
    thL[2].y = pixelData->k * (x0.z - pixelData->illuminationP.z);

    th_c = curr_constStr->aperture_C_l_plus_aperture_C_v_plus_LOG_2_PI/2 + -curr_constStr->sigt * dz;
    for (mixtureIdx = 0; mixtureIdx < curr_constStr->mixturesNum; mixtureIdx++)
    {
        gamma_s = curr_constStr->mixtureMu[mixtureIdx];
        sqrtMu = complexSqrt(
            complexSquare(cfma(gamma_s , w0.x , thL[0])) + 
            complexSquare(cfma(gamma_s , w0.y , thL[1])) +
            complexSquare(cfma(gamma_s , w0.z , thL[2])));

        log_nu = th_c + curr_constStr->mixtureC[mixtureIdx];
//         el = el + (complexExponent(sqrtMu + log_nu) - complexExponent( log_nu - sqrtMu)) / sqrtMu;
        el = el + (complexExponent(sqrtMu + log_nu)) / sqrtMu;
    }

    return el;
}

template<typename T, typename T2, typename T3>
__device__ T2 ev(T3 xb, T3 wb,
        const pixel_entry<T,T3>* __restrict__ pixelData)
{
    constStructre<T>* curr_constStr = (constStructre<T>*) &constStr_scat_mult;
    T2 ev = {0};
    T2 thV[3];
    T2 sqrtMu;
    T th_c, log_nu;

    T bd1, bd2, bd3, dz;
    unsigned char mixtureIdx;
    T gamma_s;

    // Distance to edge of sample
    bd1 = abs(pixelData->viewDir.x) < 1e-8 ? INFINITY : pixelData->viewDir.x < 0 ? (curr_constStr->box_min[0] - xb.x) / ( pixelData->viewDir.x) : (curr_constStr->box_max[0] - xb.x) / (pixelData->viewDir.x); 
    bd2 = abs(pixelData->viewDir.y) < 1e-8 ? INFINITY : pixelData->viewDir.y < 0 ? (curr_constStr->box_min[1] - xb.y) / ( pixelData->viewDir.y) : (curr_constStr->box_max[1] - xb.y) / (pixelData->viewDir.y); 
    bd3 = abs(pixelData->viewDir.z) < 1e-8 ? INFINITY : pixelData->viewDir.z < 0 ? (curr_constStr->box_min[2] - xb.z) / ( pixelData->viewDir.z) : (curr_constStr->box_max[2] - xb.z) / (pixelData->viewDir.z);
    dz  = fmin(fmin(bd1, bd2), bd3);

    // The complex throughput
    thV[0].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.x;
    thV[1].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.y;
    thV[2].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.z;

    thV[0].y = pixelData->k * (pixelData->viewP.x - xb.x);
    thV[1].y = pixelData->k * (pixelData->viewP.y - xb.y);
    thV[2].y = pixelData->k * (pixelData->viewP.z - xb.z);
         
    th_c = curr_constStr->aperture_C_l_plus_aperture_C_v_plus_LOG_2_PI/2 + ( - curr_constStr->sigt * dz);
    // For each mixture component
    for (mixtureIdx = 0; mixtureIdx < curr_constStr->mixturesNum; mixtureIdx++)
    {
        gamma_s = curr_constStr->mixtureMu[mixtureIdx];
        sqrtMu = complexSqrt(
            complexSquare(cfma(gamma_s , wb.x , thV[0])) +
            complexSquare(cfma(gamma_s , wb.y , thV[1])) +
            complexSquare(cfma(gamma_s , wb.z , thV[2])));

        log_nu = th_c + curr_constStr->mixtureC[mixtureIdx];
//         ev = ev + (complexExponent(sqrtMu + log_nu) - complexExponent(log_nu - sqrtMu)) / sqrtMu;
        ev = ev + (complexExponent(sqrtMu + log_nu)) / sqrtMu;
    }

    return ev;
}

template<typename T, typename T2, typename T3>
__device__  void scatteringLoop(ub32 current_path, ub32 current_pixel_offset,
        const T3* __restrict__ x0, const T3* __restrict__ xb,
        const T3* __restrict__ w0, const T3* __restrict__ wb,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_correlation<T,T3> *pixelData,
        const T2* __restrict__ constPath)
{
    T2 res = el<T,T2,T3>(x0[current_path],w0[current_path],&(pixelData + current_pixel_offset)->pixel_1) * 
             ev<T,T2,T3>(xb[current_path],wb[current_path],&(pixelData + current_pixel_offset)->pixel_1);
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;

    res = el<T,T2,T3>(x0[current_path],w0[current_path],&(pixelData + current_pixel_offset)->pixel_2) * 
          ev<T,T2,T3>(xb[current_path],wb[current_path],&(pixelData + current_pixel_offset)->pixel_2);
    T2 res_orig;
    res_orig.x = u_res_x[threadIdx.x]; res_orig.y = u_res_y[threadIdx.x];
    T2 u_mult = conjMult(res_orig,res);

    // Each thread copies the sum of all mixture to the shared memory
    u_mult = u_mult * constPath[current_path];
    u_res_x[threadIdx.x] = u_mult.x;
    u_res_y[threadIdx.x] = u_mult.y;
}

template<typename T, typename T2, typename T3>
__device__ void scatteringLoop(ub32 current_path, ub32 current_pixel_offset,
        const T3* __restrict__ x0, const T3* __restrict__ xb,
        const T3* __restrict__ w0, const T3* __restrict__ wb,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_field<T,T3> *pixelData,
        const T2* __restrict__ constPath)
{
    T2 res = el<T,T2,T3>(x0[current_path],w0[current_path],&(pixelData + current_pixel_offset)->pixel) *
             ev<T,T2,T3>(xb[current_path],wb[current_path],&(pixelData + current_pixel_offset)->pixel);

    res = res * constPath[current_path];
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;
}

template<typename T, typename T2, typename T3, typename correlationType, unsigned int blocksNum>
__launch_bounds__(THREADS_NUM)
__global__ void multipleScattering_kernel(T2* u,
        const T3* __restrict__ x0, const T3* __restrict__ xb,
        const T3* __restrict__ w0, const T3* __restrict__ wb,
        const correlationType* __restrict__ dataIn, const T2* __restrict__ constPath,
        const ub32* __restrict__ is_inside)
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
    ub16 current_pixel_offset = threadIdx.x % blocksNum;
    ub16 current_path_num = threadIdx.x / blocksNum;
    bool is_legal_pixel = (current_pixel_offset < pixelsPerBlock) & (current_path_num < is_inside_elemes);

    if(is_legal_pixel)
    {
        ub32 current_path = is_inside[current_path_num];
        scatteringLoop<T,T2,T3>(current_path, current_pixel_offset, x0, xb, w0, wb, u_res_x, u_res_y, pixelData, constPath);
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

template<typename T, typename T2, typename T3>
void multipleScattering(T2 *u, const T3 *x0, const T3 *xb, const T3 *w0, const T3 *wb, const T2 *constPath, const pathsList<T>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, const void* globalMem)
{
    cudaMemcpyToSymbol(is_inside_elemes, &total_paths_number, sizeof(ub32));

    if(is_correlation)
    {
        const entryStructre_correlation<T,T3> *globalMem_corr = (const entryStructre_correlation<T,T3> *) globalMem;

        if(total_paths_number > 512)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 1><<<total_elements, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 256)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 2><<<(total_elements - 1)/2 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 128)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 4><<<(total_elements - 1)/4 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 64)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 8><<<(total_elements - 1)/8 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 32)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 16><<<(total_elements - 1)/16 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 16)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 32><<<(total_elements - 1)/32 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 8)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 64><<<(total_elements - 1)/64 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
        else
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_correlation<T,T3>, 128><<<(total_elements - 1)/128 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_corr,constPath,pl->paths_numbers_list);
        }
//         mexPrintf("%d\n",cudaGetLastError());
    }
    else
    {
        const entryStructre_field<T,T3> *globalMem_field = (const entryStructre_field<T,T3> *) globalMem;

        if(total_paths_number > 512)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 1><<<total_elements, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 256)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 2><<<(total_elements - 1)/2 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 128)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 4><<<(total_elements - 1)/4 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 64)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 8><<<(total_elements - 1)/8 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 32)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 16><<<(total_elements - 1)/16 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 16)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 32><<<(total_elements - 1)/32 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else if(total_paths_number > 8)
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 64><<<(total_elements - 1)/64 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
        else
        {
            multipleScattering_kernel<T, T2, T3, entryStructre_field<T,T3>, 128><<<(total_elements - 1)/128 + 1, THREADS_NUM>>>(u,x0,xb,w0,wb,globalMem_field,constPath,pl->paths_numbers_list);
        }
    }
}

template<typename T>
void multipleScattering_setConstMem(constStructre<T> *constMem, ub32 number_of_elements)
{
    cudaMemcpyToSymbol(constStr_scat_mult, constMem, sizeof(constStructre<T>));
    cudaMemcpyToSymbol(elemNum_scat_mult, &number_of_elements, sizeof(ub32));
}

template void multipleScattering<double,double2,double3>(double2 *u, const double3 *x0, const double3 *xb, const double3 *w0, const double3 *wb, const double2 *constPath, const pathsList<double>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, const void* globalMem);

template void multipleScattering<float,float2,float3>(float2 *u, const float3 *x0, const float3 *xb, const float3 *w0, const float3 *wb, const float2 *constPath, const pathsList<float>* pl, 
        ub32 is_correlation, ub32 total_elements, ub32 total_paths_number, const void* globalMem);

template void multipleScattering_setConstMem<double>(constStructre<double> *constMem, ub32 number_of_elements);
template void multipleScattering_setConstMem<float>(constStructre<float> *constMem, ub32 number_of_elements);
