#include "nf_interface.h"
#include "nf_sampling.h"

__constant__ constStructre<double> constStr_pr_is;
__constant__ ub32 elemNum_pr_is;

template<typename T>
__global__ void preprocess_smpPositionSum_buildZ(importanceSampling<T>* IS)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_pr_is;

    ub32 z_sample_num = gridDim.x;
    ub32 z_totalSamples = blockDim.x * z_sample_num;

    T sample_w = curr_constStr->box_max[2] - curr_constStr->box_min[2];
    ub32 z_currentSample = threadIdx.x + blockIdx.x * blockDim.x;
 
    T z = curr_constStr->box_min[2] + (z_currentSample * sample_w)/((T)(z_totalSamples - 1));
    IS->z_sample[z_currentSample] = z;
}


template<typename T>
__global__ void preprocess_smpPositionSum_buildCdf(importanceSampling<T>* IS)
{
    __shared__ volatile T pdf_vals[THREADS_NUM];
    __shared__ volatile T cdf_vals[THREADS_NUM];
    __shared__ T cdf_sum_vals[THREADS_NUM];
    __shared__ T cdf_total_sum;
    __shared__ volatile T cdf_tmp_sum[32];
    // blockIdx.x: z0 grid sample
    // blockIdx.y: mixtures grid sample

    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_pr_is;

    ub32 z_sample_num = IS->z_sample_num;
    ub32 z_totalSamples = THREADS_NUM * z_sample_num;

    ub32 baseEntryNum = blockIdx.x * z_totalSamples + 
                        blockIdx.y * z_totalSamples * gridDim.x;

    // calculate pdf
    T z0 = IS->z0_sample[blockIdx.x];
    T gamma_s = curr_constStr->mixtureMu[blockIdx.y];
    T gamma_a = curr_constStr->aperture_kappa_l;

    for(ub32 zBlockNum = 0; zBlockNum < z_sample_num; zBlockNum++)
    {
        ub32 transpose_thread = (threadIdx.x % 32) * 32 + threadIdx.x / 32;
        ub32 z_currentSample = threadIdx.x + zBlockNum * blockDim.x;
        T z = IS->z_sample[z_currentSample];
        T wz_squared = (gamma_a + gamma_s) / (4 * M_PI * M_PI) + (z - z0) * (z - z0) / gamma_a;
        pdf_vals[transpose_thread] = exp(-2.0 * curr_constStr->sigt * (z - curr_constStr->box_min[2]) -log(wz_squared) );
//         pdf_vals[threadIdx.x] =exp(-2.0 * curr_constStr->sigt * (z - curr_constStr->box_min[2]));
        __syncthreads();
        // sum all elements, using wraps
        if (threadIdx.x < 32 )
        {
            cdf_vals[threadIdx.x      ] =                               pdf_vals[threadIdx.x      ];
            cdf_vals[threadIdx.x + 32 ] = cdf_vals[threadIdx.x      ] + pdf_vals[threadIdx.x + 32 ];
            cdf_vals[threadIdx.x + 64 ] = cdf_vals[threadIdx.x + 32 ] + pdf_vals[threadIdx.x + 64 ];
            cdf_vals[threadIdx.x + 96 ] = cdf_vals[threadIdx.x + 64 ] + pdf_vals[threadIdx.x + 96 ];
            cdf_vals[threadIdx.x + 128] = cdf_vals[threadIdx.x + 96 ] + pdf_vals[threadIdx.x + 128];
            cdf_vals[threadIdx.x + 160] = cdf_vals[threadIdx.x + 128] + pdf_vals[threadIdx.x + 160];
            cdf_vals[threadIdx.x + 192] = cdf_vals[threadIdx.x + 160] + pdf_vals[threadIdx.x + 192];
            cdf_vals[threadIdx.x + 224] = cdf_vals[threadIdx.x + 192] + pdf_vals[threadIdx.x + 224];
            cdf_vals[threadIdx.x + 256] = cdf_vals[threadIdx.x + 224] + pdf_vals[threadIdx.x + 256];
            cdf_vals[threadIdx.x + 288] = cdf_vals[threadIdx.x + 256] + pdf_vals[threadIdx.x + 288];
            cdf_vals[threadIdx.x + 320] = cdf_vals[threadIdx.x + 288] + pdf_vals[threadIdx.x + 320];
            cdf_vals[threadIdx.x + 352] = cdf_vals[threadIdx.x + 320] + pdf_vals[threadIdx.x + 352];
            cdf_vals[threadIdx.x + 384] = cdf_vals[threadIdx.x + 352] + pdf_vals[threadIdx.x + 384];
            cdf_vals[threadIdx.x + 416] = cdf_vals[threadIdx.x + 384] + pdf_vals[threadIdx.x + 416];
            cdf_vals[threadIdx.x + 448] = cdf_vals[threadIdx.x + 416] + pdf_vals[threadIdx.x + 448];
            cdf_vals[threadIdx.x + 480] = cdf_vals[threadIdx.x + 448] + pdf_vals[threadIdx.x + 480];
            cdf_vals[threadIdx.x + 512] = cdf_vals[threadIdx.x + 480] + pdf_vals[threadIdx.x + 512];
            cdf_vals[threadIdx.x + 544] = cdf_vals[threadIdx.x + 512] + pdf_vals[threadIdx.x + 544];
            cdf_vals[threadIdx.x + 576] = cdf_vals[threadIdx.x + 544] + pdf_vals[threadIdx.x + 576];
            cdf_vals[threadIdx.x + 608] = cdf_vals[threadIdx.x + 576] + pdf_vals[threadIdx.x + 608];
            cdf_vals[threadIdx.x + 640] = cdf_vals[threadIdx.x + 608] + pdf_vals[threadIdx.x + 640];
            cdf_vals[threadIdx.x + 672] = cdf_vals[threadIdx.x + 640] + pdf_vals[threadIdx.x + 672];
            cdf_vals[threadIdx.x + 704] = cdf_vals[threadIdx.x + 672] + pdf_vals[threadIdx.x + 704];
            cdf_vals[threadIdx.x + 736] = cdf_vals[threadIdx.x + 704] + pdf_vals[threadIdx.x + 736];
            cdf_vals[threadIdx.x + 768] = cdf_vals[threadIdx.x + 736] + pdf_vals[threadIdx.x + 768];
            cdf_vals[threadIdx.x + 800] = cdf_vals[threadIdx.x + 768] + pdf_vals[threadIdx.x + 800];
            cdf_vals[threadIdx.x + 832] = cdf_vals[threadIdx.x + 800] + pdf_vals[threadIdx.x + 832];
            cdf_vals[threadIdx.x + 864] = cdf_vals[threadIdx.x + 832] + pdf_vals[threadIdx.x + 864];
            cdf_vals[threadIdx.x + 896] = cdf_vals[threadIdx.x + 864] + pdf_vals[threadIdx.x + 896];
            cdf_vals[threadIdx.x + 928] = cdf_vals[threadIdx.x + 896] + pdf_vals[threadIdx.x + 928];
            cdf_vals[threadIdx.x + 960] = cdf_vals[threadIdx.x + 928] + pdf_vals[threadIdx.x + 960];
            cdf_vals[threadIdx.x + 992] = cdf_vals[threadIdx.x + 960] + pdf_vals[threadIdx.x + 992];

            if(threadIdx.x == 0)
            {
                if(zBlockNum == 0) {cdf_total_sum = (T) 0.0;}
                cdf_tmp_sum[0] = (T) 0.0;
                for(ub32 tmpCdfSumIter = 1; tmpCdfSumIter < 32; tmpCdfSumIter++ )
                {
                    cdf_tmp_sum[tmpCdfSumIter] = cdf_tmp_sum[tmpCdfSumIter - 1] + cdf_vals[991 + tmpCdfSumIter];
                }
                cdf_sum_vals[zBlockNum] = cdf_vals[1023] + cdf_tmp_sum[31];
                cdf_total_sum += cdf_sum_vals[zBlockNum];
            }    
        }
        __syncthreads();

        IS->sampleCdf[z_currentSample + baseEntryNum] = cdf_vals[transpose_thread] + cdf_tmp_sum[threadIdx.x / 32];
        IS->samplePdf[z_currentSample + baseEntryNum] = pdf_vals[transpose_thread];
    }

    __syncthreads();
    T addCdfVal = (T)0.0;
    T pdf_mult = ((z_totalSamples - 1)/(curr_constStr->box_max[2] - curr_constStr->box_min[2])) / cdf_total_sum;
    for(ub32 zBlockNum = 0; zBlockNum < z_sample_num; zBlockNum++)
    {
        ub32 z_currentSample = threadIdx.x + zBlockNum * blockDim.x;
        IS->samplePdf[z_currentSample + baseEntryNum] *= pdf_mult;
        IS->sampleCdf[z_currentSample + baseEntryNum] = (addCdfVal + IS->sampleCdf[z_currentSample + baseEntryNum] ) / cdf_total_sum;
        addCdfVal += cdf_sum_vals[zBlockNum];
    }
    
}

template<typename T,typename T3>
__device__ T getZ(const entryStructre_correlation<T,T3>* __restrict__ dataIn, bool isCorr)
{
    return isCorr ? dataIn->pixel_2.illuminationP.z : dataIn->pixel_1.illuminationP.z;
}

template<typename T,typename T3>
__device__ T getZ(const entryStructre_field<T,T3>* __restrict__ dataIn, bool isCorr)
{
    return dataIn->pixel.illuminationP.z;
}    


// CAN BE MADE IN LOG(N)
template<typename T, typename T3, typename correlationType>
__global__ void preprocess_smpPositionSum_buildElemantToZ0(const correlationType* __restrict__ dataIn,
        importanceSampling<T>* IS, bool isCorr)
{
    ub32 th_idx = threadIdx.x + blockIdx.x * blockDim.x;

    if(th_idx < elemNum_pr_is)
    {
        T z0 = getZ<T,T3>(dataIn + th_idx,isCorr);
        ub32 closetsIdx;
        ub32 z0_grid_num;
        bool reachToEnd = true;
        T z0_grid;

        for(z0_grid_num = 0; z0_grid_num < IS->z0_sample_num; z0_grid_num++ )
        {
            z0_grid = IS->z0_sample[z0_grid_num];
            if(z0_grid >= (z0 - 1e-8)){
                reachToEnd = false;
                break;
            }
        }

        if(z0_grid_num == 0)
        {
            closetsIdx = 0;
        }
        else if(reachToEnd)
        {
            closetsIdx = IS->z0_sample_num - 1;
        }
        else
        {
            closetsIdx = (abs(IS->z0_sample[z0_grid_num] - z0) < abs(IS->z0_sample[z0_grid_num - 1] - z0)) ?
                z0_grid_num : z0_grid_num - 1;
        }

        if(isCorr)
        {
            IS->element_to_z0_idx_2[th_idx] = closetsIdx;
        }
        else
        {
            IS->element_to_z0_idx[th_idx] = closetsIdx;
        }
    }
}

template<typename T,typename T2,typename T3>
importanceSampling<T>* sampling_preprocess(input_st<T,T2,T3> *data_in, void* gpuGlobalData)
{
    importanceSampling<T> IS_struct, *gpu_IS_struct;
    ub32 mixturesSize = data_in->scattering_vMF_mixture.mixture_size;

    if(data_in->iS.position_type == 4)
    {
        // initiate and allocate all data related to improtance sampling
        IS_struct.z0_sample_num = data_in->iS.z0_sample_num;
        IS_struct.z_sample_num = data_in->iS.z_sample_num;

        cudaMalloc(&IS_struct.samplePdf, IS_struct.z0_sample_num * IS_struct.z_sample_num * mixturesSize * THREADS_NUM * sizeof(T));
        cudaMalloc(&IS_struct.sampleCdf, IS_struct.z0_sample_num * IS_struct.z_sample_num * mixturesSize * THREADS_NUM * sizeof(T));
        cudaMalloc(&IS_struct.element_to_z0_idx, data_in->total_elements * sizeof(ub32));
        cudaMalloc(&IS_struct.z0_sample, IS_struct.z0_sample_num * sizeof(T));
        cudaMalloc(&IS_struct.z_sample, IS_struct.z_sample_num * THREADS_NUM * sizeof(T));

        if(data_in->is_correlation)
        {
            cudaMalloc(&IS_struct.element_to_z0_idx_2, data_in->total_elements * sizeof(ub32));
        }
        else
        {
            IS_struct.element_to_z0_idx_2 = 0;
        }

        T alpha_sum = 0;

        for(ub32 mixNum = 0; mixNum < mixturesSize; mixNum++)
        {
            T kappa = abs(data_in->scattering_vMF_mixture.mixture_mu[mixNum]);
            T c = data_in->scattering_vMF_mixture.mixture_c[mixNum] + log(data_in->scattering_vMF_mixture.mixture_alpha[mixNum]);
            T log_vMF_norm_factor = kappa < 1e-5 ? -log(4*M_PI) : 
                log(kappa) - (kappa > 10.0 ? log(2*M_PI) * kappa : log(4*M_PI*sinh(kappa)));
            T alpha = exp(c + log_vMF_norm_factor);
            alpha_sum = alpha_sum + alpha;
            IS_struct.alphaPdf[mixNum] = alpha;
        }

        for(ub32 mixNum = 0; mixNum < mixturesSize; mixNum++)
        {
            IS_struct.alphaPdf[mixNum] = IS_struct.alphaPdf[mixNum] / alpha_sum;
            IS_struct.alphaCdf[mixNum] = (mixNum == 0 ? (T) 0.0 : IS_struct.alphaCdf[mixNum - 1]) + IS_struct.alphaPdf[mixNum];
        }

        IS_struct.px_lines = data_in->total_elements < IS_MAX_BLOCKS ? data_in->total_elements : IS_MAX_BLOCKS;
        cudaMalloc(&IS_struct.px_table, IS_struct.px_lines * THREADS_NUM * sizeof(T));
        cudaMalloc(&IS_struct.z_idx, THREADS_NUM * sizeof(ub32));

        cudaMemcpy(IS_struct.z0_sample, data_in->iS.z0_samples, sizeof(T) * IS_struct.z0_sample_num, cudaMemcpyHostToDevice);

        // allocate the importance sampling sturctre in gpu
        cudaMalloc(&gpu_IS_struct, sizeof(importanceSampling<T>));
        cudaMemcpy(gpu_IS_struct, &IS_struct, sizeof(importanceSampling<T>), cudaMemcpyHostToDevice);

        // run preprocess procedure
        preprocess_smpPositionSum_buildZ<T><<<IS_struct.z_sample_num, THREADS_NUM>>>(gpu_IS_struct);
        preprocess_smpPositionSum_buildCdf<T><<<dim3(IS_struct.z0_sample_num,mixturesSize,1), THREADS_NUM>>>(gpu_IS_struct);

        if(data_in->is_correlation)
        {
            preprocess_smpPositionSum_buildElemantToZ0<T,T3,entryStructre_correlation<T,T3> ><<<(data_in->total_elements - 1)/THREADS_NUM + 1, THREADS_NUM>>>((entryStructre_correlation<T,T3> *)gpuGlobalData,gpu_IS_struct,false);
            preprocess_smpPositionSum_buildElemantToZ0<T,T3,entryStructre_correlation<T,T3> ><<<(data_in->total_elements - 1)/THREADS_NUM + 1, THREADS_NUM>>>((entryStructre_correlation<T,T3> *)gpuGlobalData,gpu_IS_struct,true);
        }
        else   
        {
            preprocess_smpPositionSum_buildElemantToZ0<T,T3,entryStructre_field<T,T3> ><<<(data_in->total_elements - 1)/THREADS_NUM + 1, THREADS_NUM>>>((entryStructre_field<T,T3> *)gpuGlobalData,gpu_IS_struct,false);
        }       
    }
    else
    {
        gpu_IS_struct = 0;
    }

    return gpu_IS_struct;
}

template<typename T>
void sampling_preprocess_free(importanceSampling<T> *IS_struct)
{
    if(IS_struct != 0)
    {
        importanceSampling<T>* cpu_IS = (importanceSampling<T>*) malloc(sizeof(importanceSampling<T>));
        cudaMemcpy(cpu_IS,IS_struct, sizeof(importanceSampling<T>),cudaMemcpyDeviceToHost);

        cudaFree(cpu_IS->samplePdf);
        cudaFree(cpu_IS->sampleCdf);
        cudaFree(cpu_IS->element_to_z0_idx);
        cudaFree(cpu_IS->z0_sample);
        cudaFree(cpu_IS->z_sample);
        cudaFree(cpu_IS->element_to_z0_idx_2);
        cudaFree(cpu_IS->px_table);
        cudaFree(cpu_IS->z_idx);

        cudaFree(IS_struct);
        free(cpu_IS);
    }
}

template<typename T>
void sampling_preprocess_setConstMem(constStructre<T> *constMem, ub32 number_of_elements)
{
    cudaMemcpyToSymbol(constStr_pr_is, constMem, sizeof(constStructre<T>));
    cudaMemcpyToSymbol(elemNum_pr_is, &number_of_elements, sizeof(ub32));
}

template importanceSampling<double>* sampling_preprocess<double,double2,double3>(input_st<double,double2,double3> *data_in, void* gpuGlobalData);
template importanceSampling<float>* sampling_preprocess<float,float2,float3>(input_st<float,float2,float3> *data_in, void* gpuGlobalData);

template void sampling_preprocess_free<double>(importanceSampling<double> *IS_struct);
template void sampling_preprocess_free<float>(importanceSampling<float> *IS_struct);

template void sampling_preprocess_setConstMem<double>(constStructre<double> *constMem, ub32 number_of_elements);
template void sampling_preprocess_setConstMem<float>(constStructre<float> *constMem, ub32 number_of_elements);
