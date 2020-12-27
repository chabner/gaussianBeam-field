#include "nf_header.h"
#include "nf_scattering.h"
#include "cuMath.cu"

__constant__ constStructre<double> constStr_scat_single;

template<typename T, typename T2, typename T3>
__device__ T2 singleScattering_mainLoop(const T3 currentX, const pixel_entry<T,T3>* pixelData)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_single;

    T2 thL[3];
    T2 thV[3];
    T th_c, log_nu;

    T bd1, bd2, bd3, dz;

    unsigned char mixtureIdx;

    T gamma_s, real_abs_mu;
    T2 gamma_s_over_beta_0;
    T3 w;

    T2 expVal;

    T2 C, sqrtMu;
    T2 res = {0};
        
    // l Throughput: a^{/tilde}_i

    // Distance to edge of sample
    bd1 = abs(pixelData->illuminationDir.x) < 1e-8 ? INFINITY : pixelData->illuminationDir.x >= 0 ? (currentX.x - curr_constStr->box_min[0]) / ( pixelData->illuminationDir.x) : (currentX.x - curr_constStr->box_max[0]) / (pixelData->illuminationDir.x);
    bd2 = abs(pixelData->illuminationDir.y) < 1e-8 ? INFINITY : pixelData->illuminationDir.y >= 0 ? (currentX.y - curr_constStr->box_min[1]) / ( pixelData->illuminationDir.y) : (currentX.y - curr_constStr->box_max[1]) / (pixelData->illuminationDir.y); 
    bd3 = abs(pixelData->illuminationDir.z) < 1e-8 ? INFINITY : pixelData->illuminationDir.z >= 0 ? (currentX.z - curr_constStr->box_min[2]) / ( pixelData->illuminationDir.z) : (currentX.z - curr_constStr->box_max[2]) / (pixelData->illuminationDir.z); 
        
    dz = fmin(fmin(bd1, bd2), bd3);

    // The complex throughput
    thL[0].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.x;
    thL[1].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.y;
    thL[2].x = curr_constStr->aperture_kappa_l * pixelData->illuminationDir.z;

    thL[0].y = pixelData->k * (currentX.x - pixelData->illuminationP.x);
    thL[1].y = pixelData->k * (currentX.y - pixelData->illuminationP.y);
    thL[2].y = pixelData->k * (currentX.z - pixelData->illuminationP.z);

    // v Throughput: a^{/tilde}_v

    // Distance to edge of sample
    bd1 = abs(pixelData->viewDir.x) < 1e-8 ? INFINITY : pixelData->viewDir.x < 0 ? (curr_constStr->box_min[0] - currentX.x) / ( pixelData->viewDir.x) : (curr_constStr->box_max[0] - currentX.x) / (pixelData->viewDir.x); 
    bd2 = abs(pixelData->viewDir.y) < 1e-8 ? INFINITY : pixelData->viewDir.y < 0 ? (curr_constStr->box_min[1] - currentX.y) / ( pixelData->viewDir.y) : (curr_constStr->box_max[1] - currentX.y) / (pixelData->viewDir.y); 
    bd3 = abs(pixelData->viewDir.z) < 1e-8 ? INFINITY : pixelData->viewDir.z < 0 ? (curr_constStr->box_min[2] - currentX.z) / ( pixelData->viewDir.z) : (curr_constStr->box_max[2] - currentX.z) / (pixelData->viewDir.z);
    dz += fmin(fmin(bd1, bd2), bd3);

    // The complex throughput
    thV[0].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.x;
    thV[1].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.y;
    thV[2].x = curr_constStr->aperture_kappa_v * pixelData->viewDir.z;

    thV[0].y = pixelData->k * (pixelData->viewP.x - currentX.x);
    thV[1].y = pixelData->k * (pixelData->viewP.y - currentX.y);
    thV[2].y = pixelData->k * (pixelData->viewP.z - currentX.z);
 
    th_c = curr_constStr->aperture_C_l_plus_aperture_C_v_plus_LOG_2_PI - curr_constStr->sigt * dz;
        
    // For each mixture component
    for (mixtureIdx = 0; mixtureIdx < curr_constStr->mixturesNum; mixtureIdx++)
    {
        // Convolution with the illumination throughput
        gamma_s = curr_constStr->mixtureMu[mixtureIdx];
        if(abs(gamma_s) < 0.000000001)
        {
            gamma_s = 0.000000001;
        }
        
        gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
                complexSquare(cfma(gamma_s , pixelData->viewDir.x , thL[0])) +
                complexSquare(cfma(gamma_s , pixelData->viewDir.y , thL[1])) +
                complexSquare(cfma(gamma_s , pixelData->viewDir.z , thL[2])));

        sqrtMu = complexSqrt(
            complexSquare(cfma(gamma_s_over_beta_0 , thL[0] , thV[0])) +
            complexSquare(cfma(gamma_s_over_beta_0 , thL[1] , thV[1])) +
            complexSquare(cfma(gamma_s_over_beta_0 , thL[2] , thV[2])));

        real_abs_mu = rnorm3d(realMult(gamma_s_over_beta_0 , thL[0]),realMult(gamma_s_over_beta_0 , thL[1]),realMult(gamma_s_over_beta_0 , thL[2]));
        w.x = realMult(gamma_s_over_beta_0 , thL[0]); w.y = realMult(gamma_s_over_beta_0 , thL[1]); w.z = realMult(gamma_s_over_beta_0 , thL[2]);
//         expVal = expVal - real_abs_mu * gamma_s_over_beta_0 * cfma(w.x , thL[0] , cfma(w.y , thL[1] , w.z * thL[2]));
        gamma_s *= real_abs_mu;

        C = complexSqrt(
            complexSquare(cfma(gamma_s , w.x , thL[0])) +
            complexSquare(cfma(gamma_s , w.y , thL[1])) +
            complexSquare(cfma(gamma_s , w.z , thL[2])));

        expVal = sqrtMu + C - real_abs_mu * gamma_s_over_beta_0 * cfma(w.x , thL[0] , cfma(w.y , thL[1] , w.z * thL[2]));
        log_nu = th_c + curr_constStr->mixtureC[mixtureIdx];
                
        // integrate
//         res = res + (complexExponent(expVal + log_nu) - complexExponent(log_nu - expVal)) / (C * sqrtMu);
        res = res + (complexExponent(expVal + log_nu)) / (C * sqrtMu);
    }

    return res;
}
        
template<typename T, typename T2, typename T3>
__device__  void singleScatteringLoop(const T3* __restrict__ x0,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_correlation<T,T3> *pixelData,
        const T2* __restrict__ constPath)
{
    T2 res = singleScattering_mainLoop<T,T2,T3>(x0[threadIdx.x],&pixelData->pixel_1);
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;

    res = singleScattering_mainLoop<T,T2,T3>(x0[threadIdx.x],&pixelData->pixel_2);
    T2 res_orig;
    res_orig.x = u_res_x[threadIdx.x]; res_orig.y = u_res_y[threadIdx.x];
    T2 u_mult = conjMult(res_orig,res);

    // Each thread copies the sum of all mixture to the shared memory
    u_mult = u_mult * constPath[threadIdx.x];
    u_res_x[threadIdx.x] = u_mult.x;
    u_res_y[threadIdx.x] = u_mult.y;
}

template<typename T, typename T2, typename T3>
__device__  void singleScatteringLoop(const T3* __restrict__ x0,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_field<T,T3> *pixelData,
        const T2* __restrict__ constPath)
{
    T2 res = singleScattering_mainLoop<T,T2,T3>(x0[threadIdx.x],&pixelData->pixel);
    res = res * constPath[threadIdx.x];
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;
}
        
template<typename T, typename T2, typename T3, typename correlationType>
__launch_bounds__(THREADS_NUM)
__global__ void singleScattering_kernel(T2* u, const T3* __restrict__ x0,
        const correlationType* __restrict__ dataIn, const T2* __restrict__ constPath)
{
    __shared__ volatile T u_res_x[THREADS_NUM];
    __shared__ volatile T u_res_y[THREADS_NUM];
    __shared__ correlationType pixelData;

    // copy the relevant data to shared memory
    int* iDest = (int*)&pixelData;
    const int* iSrc = (const int*)(dataIn + blockIdx.x);

    if(threadIdx.x < sizeof(correlationType) / sizeof(int))
    {
        iDest[threadIdx.x] = iSrc[threadIdx.x];
    }

    __syncthreads();
    singleScatteringLoop<T,T2,T3>(x0,u_res_x,u_res_y,&pixelData,constPath);
    __syncthreads();

    // Reduce sum of all shared memory to a single result
    if (threadIdx.x < 512){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 512];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 512];
    }  __syncthreads();
    if (threadIdx.x < 256){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 256];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 256];
    }  __syncthreads();
    if (threadIdx.x < 128){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 128];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 128];
    }  __syncthreads();
    if (threadIdx.x < 64 ){
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 64 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 64 ];
    }  __syncthreads();

    if (threadIdx.x < 32 ) // warpReduce
    {
        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 32];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 32];

        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 16];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 16];

        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 8 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 8 ];

        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 4 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 4 ];

        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 2 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 2 ];

        u_res_x[threadIdx.x] = u_res_x[threadIdx.x] + u_res_x[threadIdx.x + 1 ];
        u_res_y[threadIdx.x] = u_res_y[threadIdx.x] + u_res_y[threadIdx.x + 1 ];
    }

//     copy the final result to the global mem
    if (threadIdx.x == 0) { u[blockIdx.x].x += u_res_x[0]; u[blockIdx.x].y += u_res_y[0]; }
}

template<typename T, typename T2, typename T3>
void singleScattering(T2 *us, const T3 *x0, const T2* constPath, ub32 is_correlation, ub32 total_elements, const void* globalMem)
{
    if(is_correlation)
    {
        const entryStructre_correlation<T,T3> *globalMem_corr = (const entryStructre_correlation<T,T3> *) globalMem;
        singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3> ><<<total_elements, THREADS_NUM>>>(us,x0,globalMem_corr,constPath);
    }
    else
    {
        const entryStructre_field<T,T3> *globalMem_field = (const entryStructre_field<T,T3> *) globalMem;
        singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3> ><<<total_elements, THREADS_NUM>>>(us,x0,globalMem_field,constPath);
    }
}

template<typename T>
void singleScattering_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_scat_single, constMem, sizeof(constStructre<T>));
}

template void singleScattering<double,double2,double3>(double2 *us, const double3 *x0, const double2* constPath, ub32 is_correlation, ub32 total_elements, const void* globalMem);
template void singleScattering<float,float2,float3>(float2 *us, const float3 *x0, const float2* constPath, ub32 is_correlation, ub32 total_elements, const void* globalMem);

template void singleScattering_setConstMem<double>(constStructre<double> *constMem);
template void singleScattering_setConstMem<float>(constStructre<float> *constMem);
