#include "nf_header.h"
#include "nf_sampling.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ constStructre<double> constStr_smp_pos;
__constant__ ub32 elemNum_smp_pos;
__constant__ ub32 pixel_shift;

template<typename T, typename T3>
__device__ pixel_entry<T,T3> getPxixelEntry(const entryStructre_correlation<T,T3>* __restrict__ dataIn, ub32 n, ub32 c)
{
    return c == 0 ? dataIn[n].pixel_1 : dataIn[n].pixel_2;
}

template<typename T, typename T3>
__device__ pixel_entry<T,T3> getPxixelEntry(const entryStructre_field<T,T3>* __restrict__ dataIn, ub32 n, ub32 c)
{
    return dataIn[n].pixel;
}

template<typename T>
__device__ T calcP_pixel(const importanceSampling<T>* IS, T dx_sq_plus_dy_sq, T z, ub32 z0_idx, ub32 z_idx)
{    
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    ub32 totalSamples = IS->z_sample_num * THREADS_NUM;

    T z0 = IS->z0_sample[z0_idx];
    T p0 = 0;
    T gamma_a = curr_constStr->aperture_kappa_l;

    for(ub32 current_mixtureNum = 0; current_mixtureNum < curr_constStr->mixturesNum; current_mixtureNum++)
    {
        ub32 sampleSection = totalSamples * (z0_idx + IS->z0_sample_num * current_mixtureNum);

        T gamma_s = curr_constStr->mixtureMu[current_mixtureNum];
        T wz_squared = (gamma_a + gamma_s) / (4 * M_PI * M_PI) + (z - z0) * (z - z0) / gamma_a;

        T pz = IS->samplePdf[sampleSection + z_idx];
        T px_py = exp(-(dx_sq_plus_dy_sq) / (2 * wz_squared) ) / ((2 * wz_squared) * M_PI);
        p0 += pz * px_py * IS->alphaPdf[current_mixtureNum];
    }              

    return p0;
}


template<typename T, typename T3>
__device__ T calcP(const entryStructre_correlation<T,T3>* __restrict__ dataIn,
        const importanceSampling<T>* IS, T3 curr_x)
{
    // in correlation we calculate the mean of p1 and p2
            
    ub32 current_pixel = pixel_shift + blockIdx.x;

    T dx   = curr_x.x - (dataIn + current_pixel)->pixel_1.illuminationP.x;
    T dy   = curr_x.y - (dataIn + current_pixel)->pixel_1.illuminationP.y;
    T dx_2 = curr_x.x - (dataIn + current_pixel)->pixel_2.illuminationP.x;
    T dy_2 = curr_x.y - (dataIn + current_pixel)->pixel_2.illuminationP.y;

    T dx_sq_plus_dy_sq   = dx   * dx   + dy   * dy;
    T dx_sq_plus_dy_sq_2 = dx_2 * dx_2 + dy_2 * dy_2;

    ub32 z0_idx   = IS->element_to_z0_idx[current_pixel];
    ub32 z0_idx_2 = IS->element_to_z0_idx[current_pixel];

    ub32 z_idx = IS->z_idx[threadIdx.x];

    return 0.5 * (calcP_pixel<T>(IS,dx_sq_plus_dy_sq  ,curr_x.z,z0_idx  ,z_idx) + 
                  calcP_pixel<T>(IS,dx_sq_plus_dy_sq_2,curr_x.z,z0_idx_2,z_idx));
}

template<typename T, typename T3>
__device__ T calcP(const entryStructre_field<T,T3>* __restrict__ dataIn,
        const importanceSampling<T>* IS, T3 curr_x)
{
    ub32 current_pixel = pixel_shift + blockIdx.x;

    T dx = curr_x.x - (dataIn + current_pixel)->pixel.illuminationP.x;
    T dy = curr_x.y - (dataIn + current_pixel)->pixel.illuminationP.y;
    T dx_sq_plus_dy_sq = dx * dx + dy * dy;
    ub32 z0_idx = IS->element_to_z0_idx[current_pixel];

    ub32 z_idx = IS->z_idx[threadIdx.x];

    return calcP_pixel<T>(IS,dx_sq_plus_dy_sq,curr_x.z,z0_idx,z_idx);
}

template<typename T, typename T3, typename correlationType, bool isCorr>
__launch_bounds__(32)
__global__ void samplePosition_gbSum_sample(ub32* x_rep, T3* x0, ub32* n, ub32* c,
        const correlationType* __restrict__ dataIn, importanceSampling<T>* IS, time_t seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    __shared__ ub32 local_x_rep;

    if(threadIdx.x == 0)
    {
        local_x_rep = 0;
    }
        
    __syncthreads();

    ub32 th_id = threadIdx.x + (blockDim.x * blockIdx.x);

    curandState_t state;
    curand_init(seed, th_id, 0, &state);

    ub32 current_n = curand(&state) % elemNum_smp_pos;
    ub32 current_c = (isCorr ? curand(&state) % 2 : 0);

    // randomize alpha
    T alphaRand;
    ub32 current_alpha_idx;

    // randomize sample
    ub32 totalSamples = IS->z_sample_num * THREADS_NUM;
    T x, y, z;
    ub32 m;
    bool is_inside;
    ub32 z0_sample_num;
    pixel_entry<T,T3> currPixel;

    do{
        current_alpha_idx = 0;
        randomUniform(&state,&alphaRand);
        while(IS->alphaCdf[current_alpha_idx] < alphaRand && current_alpha_idx < curr_constStr->mixturesNum)
        {
            current_alpha_idx++;
        }

        // z sample
        if(current_c == 0)
        {
            z0_sample_num = IS->element_to_z0_idx[current_n];
        }
        else
        {
            z0_sample_num = IS->element_to_z0_idx_2[current_n];
        }

        ub32 sampleSection = totalSamples *                     z0_sample_num + 
                                      totalSamples * IS->z0_sample_num * current_alpha_idx;

        T zRand;
        randomUniform(&state,&zRand);

        m = binary_search<T>(IS->sampleCdf + sampleSection, zRand, totalSamples);
        z = IS->z_sample[m];

        // x-y sample

        currPixel = getPxixelEntry<T,T3>(dataIn, current_n, current_c);

        T z0 = IS->z0_sample[z0_sample_num];

        // not sure if abs is right here for negative gamma_s
        T gamma_s = abs(curr_constStr->mixtureMu[current_alpha_idx]);
        T gamma_a = curr_constStr->aperture_kappa_l;
        T wz = sqrt((gamma_a + gamma_s) / (4 * M_PI * M_PI) + (z - z0) * (z - z0) / gamma_a);

        T rand_x, rand_y;
        randomNormal(&state, &rand_x);
        randomNormal(&state, &rand_y);
        x = wz * rand_x + currPixel.illuminationP.x;
        y = wz * rand_y + currPixel.illuminationP.y;

        is_inside = (x >= curr_constStr->box_min[0] && x <= curr_constStr->box_max[0]) &&
                    (y >= curr_constStr->box_min[1] && y <= curr_constStr->box_max[1]) &&
                    (z >= curr_constStr->box_min[2] && z <= curr_constStr->box_max[2]);

        if(!is_inside)
        {
            #if (__CUDA_ARCH__ >= 600)
                atomicAdd_block(&local_x_rep,1);
            #else
                atomicAdd(&local_x_rep,1);
            #endif
            current_n = curand(&state) % elemNum_smp_pos;
            current_c = (isCorr ? curand(&state) % 2 : 0);
        }

    }
    while(!is_inside);

    x0[th_id].x = x;
    x0[th_id].y = y;
    x0[th_id].z = z;

    if(threadIdx.x == 0) {atomicAdd(x_rep,local_x_rep);} 

    n[th_id] = current_n;
    c[th_id] = current_c;

    IS->z_idx[th_id] = m;
}

template<typename T, typename T3, typename correlationType>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_gbSum_probability(const T3* __restrict__ x,
        const correlationType* __restrict__ dataIn, importanceSampling<T>* IS)
{
    IS->px_table[threadIdx.x + (blockDim.x * blockIdx.x)] = calcP<T,T3>(dataIn,IS,x[threadIdx.x]);
}

template<typename T>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_gbSum_sum(T* px0, const importanceSampling<T>* IS, ub32 callBlocks)
{
    T px = (T) 0.0;
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;

    for(ub32 iter = 0; iter < callBlocks; iter++)
    {
        px += IS->px_table[threadIdx.x + (blockDim.x * iter)];
    }
    px0[threadIdx.x] += (px  * curr_constStr->V) / elemNum_smp_pos;
}

template<typename T, typename T3, bool isCorr>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_uniform(T3* gpu_x, T* gpu_x0p, ub32* n, ub32* c, ub64 seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    T x, y, z;

    randomUniform(&state,&x);
    randomUniform(&state,&y);
    randomUniform(&state,&z);

    gpu_x[threadIdx.x].x = fma(x,curr_constStr->box_max[0] - curr_constStr->box_min[0],curr_constStr->box_min[0]);
    gpu_x[threadIdx.x].y = fma(y,curr_constStr->box_max[1] - curr_constStr->box_min[1],curr_constStr->box_min[1]);
    gpu_x[threadIdx.x].z = fma(z,curr_constStr->box_max[2] - curr_constStr->box_min[2],curr_constStr->box_min[2]);
    gpu_x0p[threadIdx.x] = (T) 1.0; 

    ub32 current_n = curand(&state) % elemNum_smp_pos;
    ub32 current_c = (isCorr ? curand(&state) % 2 : 0);

    n[threadIdx.x] = current_n;
    c[threadIdx.x] = current_c;
}

template<typename T, typename T3, bool isCorr>
__launch_bounds__(THREADS_NUM)
__global__ void samplePosition_exponential(T3* gpu_x, T* gpu_x0p, ub32* n, ub32* c, time_t seed)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_pos;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    // z is distributed exponentially
    T sigt = curr_constStr->sigt * ((T)2.0);
    T dmax = curr_constStr->box_max[2] - curr_constStr->box_min[2];
    T rmax = 1-exp(-dmax*sigt);

    T rand_z;
    randomUniform(&state,&rand_z);
    T d = -log(fma(-rand_z,rmax,(T)1.0))/(sigt);

    gpu_x[threadIdx.x].z = d + curr_constStr->box_min[2];
    gpu_x0p[threadIdx.x] = (exp(-sigt * d) / (rmax)) * sigt * dmax;

    T x, y;

    randomUniform(&state,&x);
    randomUniform(&state,&y);

    // x and y are distributed uniformly
    gpu_x[threadIdx.x].x = fma(x,curr_constStr->box_max[0] - curr_constStr->box_min[0],curr_constStr->box_min[0]);
    gpu_x[threadIdx.x].y = fma(y,curr_constStr->box_max[1] - curr_constStr->box_min[1],curr_constStr->box_min[1]);

    ub32 current_n = curand(&state) % elemNum_smp_pos;
    ub32 current_c = (isCorr ? curand(&state) % 2 : 0);

    n[threadIdx.x] = current_n;
    c[threadIdx.x] = current_c;

}

template<typename T, typename T3>
void samplePosition(T3 *x0, T* px0, ub32 *n, ub32 *c, ub32 *x_rep,
        ub32 sample_position_flag, ub32 is_correlation, ub32 total_elements,
        importanceSampling<T>* IS_struct, const void* globalMem, ub64 curr_seed, T min_prob)
{
    if(sample_position_flag == 1)
    {        
        // unifrom sampling
        if(is_correlation)
        {
            samplePosition_uniform<T,T3,true><<<1, THREADS_NUM>>>(x0,px0,n,c,curr_seed);
        }
        else
        {
            samplePosition_uniform<T,T3,false><<<1, THREADS_NUM>>>(x0,px0,n,c,curr_seed);
        }

    }
    else if(sample_position_flag == 2)
    {
        // exponential sampling
        if(is_correlation)
        {
            samplePosition_exponential<T,T3,true><<<1, THREADS_NUM>>>(x0,px0,n,c,curr_seed);
        }
        else
        {
            samplePosition_exponential<T,T3,false><<<1, THREADS_NUM>>>(x0,px0,n,c,curr_seed);

            // in field we take the sqrt of px0
            global_sqrt<T><<<1, THREADS_NUM>>>(px0);
        }
    }
    else if(sample_position_flag == 4)
    {
        // Gaussian beam sum samplimng
        cudaMemset(px0, 0, THREADS_NUM * sizeof(T));

        if(is_correlation)
        {
            entryStructre_correlation<T,T3> *globalMem_corr = (entryStructre_correlation<T,T3> *) globalMem;

            // First let's sample the points
            samplePosition_gbSum_sample<T,T3,entryStructre_correlation<T,T3>,true><<<32, THREADS_NUM/32>>>(
                    x_rep,x0,n,c,globalMem_corr,IS_struct,curr_seed);

            // Now calculate the probability of the points
            for(ub32 px_iter = 0; px_iter < total_elements; px_iter += IS_MAX_BLOCKS)
            {
                cudaMemcpyToSymbol(pixel_shift, &px_iter, sizeof(ub32));
                
                ub32 callBlocks = total_elements - px_iter > IS_MAX_BLOCKS ? IS_MAX_BLOCKS : (total_elements - px_iter);

                samplePosition_gbSum_probability<T,T3,entryStructre_correlation<T,T3>><<<callBlocks, THREADS_NUM>>>(
                        x0,globalMem_corr,IS_struct);

                samplePosition_gbSum_sum<T><<<1, THREADS_NUM>>>(px0,IS_struct,callBlocks);
            }
        }
        else
        {
            entryStructre_field<T,T3> *globalMem_field = (entryStructre_field<T,T3> *) globalMem;

            // First let's sample the points
            samplePosition_gbSum_sample<T,T3,entryStructre_field<T,T3>,false><<<32, THREADS_NUM/32>>>(
                    x_rep,x0,n,c,globalMem_field,IS_struct,curr_seed);

            // Now calculate the probability of the points
            for(ub32 px_iter = 0; px_iter < total_elements; px_iter += IS_MAX_BLOCKS)
            {
                cudaMemcpyToSymbol(pixel_shift, &px_iter, sizeof(ub32));
                
                ub32 callBlocks = total_elements - px_iter > IS_MAX_BLOCKS ? IS_MAX_BLOCKS : (total_elements - px_iter);

                samplePosition_gbSum_probability<T,T3,entryStructre_field<T,T3>><<<callBlocks, THREADS_NUM>>>(
                        x0,globalMem_field,IS_struct);

                samplePosition_gbSum_sum<T><<<1, THREADS_NUM>>>(px0,IS_struct,callBlocks);
            }

            // in field we take the sqrt of px0
            global_sqrt<T><<<1, THREADS_NUM>>>(px0);
        }
    }

    global_min<T><<<1, THREADS_NUM>>>(px0, min_prob);
}

template<typename T>
void samplePosition_setConstMem(constStructre<T> *constMem, ub32 number_of_elements)
{
    cudaMemcpyToSymbol(constStr_smp_pos, constMem, sizeof(constStructre<T>));
    cudaMemcpyToSymbol(elemNum_smp_pos, &number_of_elements, sizeof(ub32));
}

template void samplePosition<double,double3>(double3 *x0, double* px0, ub32 *n, ub32 *c, ub32 *x_rep,
        ub32 sample_position_flag, ub32 is_correlation, ub32 total_elements,
        importanceSampling<double>* IS_struct, const void* globalMem, ub64 curr_seed, double min_prob);

template void samplePosition<float,float3>(float3 *x0, float* px0, ub32 *n, ub32 *c, ub32 *x_rep,
        ub32 sample_position_flag, ub32 is_correlation, ub32 total_elements,
        importanceSampling<float>* IS_struct, const void* globalMem, ub64 curr_seed, float min_prob);

template void samplePosition_setConstMem<double>(constStructre<double> *constMem, ub32 number_of_elements);
template void samplePosition_setConstMem<float>(constStructre<float> *constMem, ub32 number_of_elements);
