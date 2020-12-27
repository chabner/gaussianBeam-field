#include "ff_header.h"
#include "ff_scattering.h"
#include "cuMath.cu"

__constant__ constStructre<double> constStr_scat_single;

template<typename T, typename T3, typename scatteringType, ub32 dims>
__device__  void calc_scatter_contribution(T *scatContrb,
        const scatteringType* scattering_struct,
        const entryStructre_correlation<T,T3>* entry_struct)
{
    scatContrb[0] = scattering_contribution<T,T3,dims,true>(scattering_struct, 
            &entry_struct->pixel_1.illumination_direction, &entry_struct->pixel_1.view_direction);

    scatContrb[1] = scattering_contribution<T,T3,dims,true>(scattering_struct, 
            &entry_struct->pixel_2.illumination_direction, &entry_struct->pixel_2.view_direction);
}

template<typename T, typename T3, typename scatteringType, ub32 dims>
__device__  void calc_scatter_contribution(T *scatContrb,
        const scatteringType* scattering_struct,
        const entryStructre_field<T,T3>* entry_struct)
{
    *scatContrb = scattering_contribution<T,T3,dims,true>(scattering_struct, 
            &entry_struct->pixel.illumination_direction, &entry_struct->pixel.view_direction);
}

template<typename T, typename T2, typename T3, ub32 dims>
__device__ T2 singleScattering_mainLoop(const T3 currentX0, const pixel_entry<T,T3>* pixelData)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_scat_single;

    T2 e;

    T bd1, bd2, dz;

    bd1 = abs(pixelData->illumination_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.x >= 0 ? (currentX0.x - curr_constStr->box_min[0]) / ( pixelData->illumination_attenuation_direction.x) : (currentX0.x - curr_constStr->box_max[0]) / (pixelData->illumination_attenuation_direction.x);
    bd2 = abs(pixelData->illumination_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.y >= 0 ? (currentX0.y - curr_constStr->box_min[1]) / ( pixelData->illumination_attenuation_direction.y) : (currentX0.y - curr_constStr->box_max[1]) / (pixelData->illumination_attenuation_direction.y); 

    if(dims == 3)
    {
        T bd3 = abs(pixelData->illumination_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->illumination_attenuation_direction.z >= 0 ? (currentX0.z - curr_constStr->box_min[2]) / ( pixelData->illumination_attenuation_direction.z) : (currentX0.z - curr_constStr->box_max[2]) / (pixelData->illumination_attenuation_direction.z);
        dz = fmin(fmin(bd1, bd2), bd3);
    }
    if(dims == 2)
    {
        dz = fmin(bd1, bd2);
    }

    bd1 = abs(pixelData->view_attenuation_direction.x) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.x < 0 ? (curr_constStr->box_min[0] - currentX0.x) / ( pixelData->view_attenuation_direction.x) : (curr_constStr->box_max[0] - currentX0.x) / (pixelData->view_attenuation_direction.x); 
    bd2 = abs(pixelData->view_attenuation_direction.y) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.y < 0 ? (curr_constStr->box_min[1] - currentX0.y) / ( pixelData->view_attenuation_direction.y) : (curr_constStr->box_max[1] - currentX0.y) / (pixelData->view_attenuation_direction.y); 

    if(dims == 3)
    {
        T bd3 = abs(pixelData->view_attenuation_direction.z) < 1e-8 ? INFINITY : pixelData->view_attenuation_direction.z < 0 ? (curr_constStr->box_min[2] - currentX0.z) / ( pixelData->view_attenuation_direction.z) : (curr_constStr->box_max[2] - currentX0.z) / (pixelData->view_attenuation_direction.z);
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
                  (currentX0.x * (pixelData->illumination_direction.x - pixelData->view_direction.x) +
                   currentX0.y * (pixelData->illumination_direction.y - pixelData->view_direction.y) +
                   currentX0.z * (pixelData->illumination_direction.z - pixelData->view_direction.z) );
    }
    if(dims == 2)
    {
        e.y = pixelData->k * 
                  (currentX0.x * (pixelData->illumination_direction.x - pixelData->view_direction.x) +
                   currentX0.y * (pixelData->illumination_direction.y - pixelData->view_direction.y) );
    }

    return complexExponent(e);
}
        
template<typename T, typename T2, typename T3, ub32 dims>
__device__  void singleScatteringLoop(const T3* x0,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_correlation<T,T3> *pixelData,
        const T2* constPath, const T* scatContrb)
{
    T2 res_1 = singleScattering_mainLoop<T,T2,T3,dims>(x0[threadIdx.x],&pixelData->pixel_1);
    T2 res_2 = singleScattering_mainLoop<T,T2,T3,dims>(x0[threadIdx.x],&pixelData->pixel_2);
    T2 u_mult = conjMult(scatContrb[0] * res_1,scatContrb[1] * res_2);

    // Each thread copies the sum of all mixture to the shared memory
    u_mult = u_mult * constPath[threadIdx.x];
    u_res_x[threadIdx.x] = u_mult.x;
    u_res_y[threadIdx.x] = u_mult.y;
}

template<typename T, typename T2, typename T3, ub32 dims>
__device__  void singleScatteringLoop(const T3* x0,
        volatile T *u_res_x, volatile T *u_res_y, entryStructre_field<T,T3> *pixelData,
        const T2* constPath, const T* scatContrb)
{
    T2 res = singleScattering_mainLoop<T,T2,T3,dims>(x0[threadIdx.x],&pixelData->pixel);
    res = scatContrb[0] * res * constPath[threadIdx.x];
    u_res_x[threadIdx.x] = res.x;
    u_res_y[threadIdx.x] = res.y;
}
        
template<typename T, typename T2, typename T3, typename correlationType, typename scatteringType, ub32 dims>
__launch_bounds__(THREADS_NUM)
__global__ void singleScattering_kernel(T2* us, const T3* x0,
        const correlationType* dataIn, const T2* constPath,
        const scatteringType* scatteringStruct)
{
    __shared__ volatile T u_res_x[THREADS_NUM];
    __shared__ volatile T u_res_y[THREADS_NUM];
    __shared__ correlationType pixelData;
    __shared__ T scatContrb[2];

    // copy the relevant data to shared memory
    int* iDest = (int*)&pixelData;
    const int* iSrc = (const int*)(dataIn + blockIdx.x);

    if(threadIdx.x < sizeof(correlationType) / sizeof(int))
    {
        iDest[threadIdx.x] = iSrc[threadIdx.x];
    }

    __syncthreads();

    if(threadIdx.x == 0)
    {
        calc_scatter_contribution<T,T3,scatteringType,dims>(scatContrb,scatteringStruct,&pixelData);
    }

    __syncthreads();
    singleScatteringLoop<T,T2,T3,dims>(x0,u_res_x,u_res_y,&pixelData,constPath,scatContrb);
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
    if (threadIdx.x == 0) { us[blockIdx.x].x += u_res_x[0]; us[blockIdx.x].y += u_res_y[0]; }
}

template<typename T, typename T2, typename T3>
void singleScattering(T2 *us, const T3 *x0, const T2* constPath, ub32 is_correlation, ub32 total_elements,
        ub32 dims, const void* globalMem, ub32 scattering_paramenter, const void* scatteringStruct)
{    
    if(is_correlation)
    {   
        if(scattering_paramenter == 1)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(const entryStructre_correlation<T,T3>*) globalMem,constPath,(const randomScatter*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,randomScatter,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(const entryStructre_correlation<T,T3>*) globalMem,constPath,(const randomScatter*) scatteringStruct);
            }
        }

        if(scattering_paramenter == 2)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_correlation<T,T3>*) globalMem,constPath,(tabulatedScatter<T>*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,tabulatedScatter<T>,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_correlation<T,T3>*) globalMem,constPath,(tabulatedScatter<T>*) scatteringStruct);
            }
        }

        if(scattering_paramenter == 3)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_correlation<T,T3>*) globalMem,constPath,(HGscatter<T>*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_correlation<T,T3>,HGscatter<T>,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_correlation<T,T3>*) globalMem,constPath,(HGscatter<T>*) scatteringStruct);
            }
        }
    }
    else
    {
        if(scattering_paramenter == 1)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,randomScatter,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(const entryStructre_field<T,T3>*) globalMem,constPath,(const randomScatter*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,randomScatter,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(const entryStructre_field<T,T3>*) globalMem,constPath,(const randomScatter*) scatteringStruct);
            }
        }

        if(scattering_paramenter == 2)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_field<T,T3>*) globalMem,constPath,(tabulatedScatter<T>*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,tabulatedScatter<T>,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_field<T,T3>*) globalMem,constPath,(tabulatedScatter<T>*) scatteringStruct);
            }
        }

        if(scattering_paramenter == 3)
        {
            if(dims == 2)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,2><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_field<T,T3>*) globalMem,constPath,(HGscatter<T>*) scatteringStruct);
            }
            if(dims == 3)
            {
                singleScattering_kernel<T,T2,T3,entryStructre_field<T,T3>,HGscatter<T>,3><<<total_elements, THREADS_NUM>>>
                        (us,x0,(entryStructre_field<T,T3>*) globalMem,constPath,(HGscatter<T>*) scatteringStruct);
            }
        }
    }
}

template<typename T>
void singleScattering_setConstMem(constStructre<T> *constMem)
{
    cudaMemcpyToSymbol(constStr_scat_single, constMem, sizeof(constStructre<T>));
}

template void singleScattering_setConstMem<double>(constStructre<double> *constMem);
template void singleScattering_setConstMem<float>(constStructre<float> *constMem);

template void singleScattering<double,double2,double3>(double2 *us, const double3 *x0, const double2* constPath, ub32 is_correlation, ub32 total_elements,
        ub32 dims, const void* globalMem, ub32 scattering_paramenter, const void* scatteringStruct);
template void singleScattering<float,float2,float3>(float2 *us, const float3 *x0, const float2* constPath, ub32 is_correlation, ub32 total_elements,
        ub32 dims, const void* globalMem, ub32 scattering_paramenter, const void* scatteringStruct);
