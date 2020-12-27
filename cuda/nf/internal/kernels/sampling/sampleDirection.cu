#include "nf_header.h"
#include "nf_sampling.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ constStructre<double> constStr_smp_dir;
__constant__ ub32 elemNum_smp_dir;
__constant__ bool isPreRand;

template<typename T, typename T2, typename T3>
__device__ void buildTh(T2 *thL, T3 *dirv, curandState_t* state,
        const T3* __restrict__ x, const ub32* __restrict__ n, const ub32* __restrict__ c,
        const entryStructre_correlation<T,T3>* __restrict__ dataIn)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_dir;
    ub32 th_id = threadIdx.x + (blockDim.x * blockIdx.x);
    T3 currentX = x[th_id];
    ub32 current_n, current_c;
            
    if(isPreRand)
    {
        current_n = n[th_id];
        current_c = c[th_id];
    }
    else
    {
        current_n = curand(state) % elemNum_smp_dir;
        current_c = curand(state) % 2;
    }

    const entryStructre_correlation<T,T3>* pixelData = dataIn + current_n;

    if(current_c == 0)
    {
        thL[0].x = curr_constStr->aperture_kappa_l * pixelData->pixel_1.illuminationDir.x;
        thL[1].x = curr_constStr->aperture_kappa_l * pixelData->pixel_1.illuminationDir.y;
        thL[2].x = curr_constStr->aperture_kappa_l * pixelData->pixel_1.illuminationDir.z;

        thL[0].y = pixelData->pixel_1.k * (currentX.x - pixelData->pixel_1.illuminationP.x);
        thL[1].y = pixelData->pixel_1.k * (currentX.y - pixelData->pixel_1.illuminationP.y);
        thL[2].y = pixelData->pixel_1.k * (currentX.z - pixelData->pixel_1.illuminationP.z);

        dirv->x = pixelData->pixel_1.viewDir.x;
        dirv->y = pixelData->pixel_1.viewDir.y;
        dirv->z = pixelData->pixel_1.viewDir.z;
    }
    else
    {
        thL[0].x = curr_constStr->aperture_kappa_l * pixelData->pixel_2.illuminationDir.x;
        thL[1].x = curr_constStr->aperture_kappa_l * pixelData->pixel_2.illuminationDir.y;
        thL[2].x = curr_constStr->aperture_kappa_l * pixelData->pixel_2.illuminationDir.z;

        thL[0].y = pixelData->pixel_2.k * (currentX.x - pixelData->pixel_2.illuminationP.x);
        thL[1].y = pixelData->pixel_2.k * (currentX.y - pixelData->pixel_2.illuminationP.y);
        thL[2].y = pixelData->pixel_2.k * (currentX.z - pixelData->pixel_2.illuminationP.z);

        dirv->x = pixelData->pixel_2.viewDir.x;
        dirv->y = pixelData->pixel_2.viewDir.y;
        dirv->z = pixelData->pixel_2.viewDir.z;
    }
}

template<typename T, typename T2, typename T3>
__device__ void buildTh(T2 *thL, T3 *dirv, curandState_t* state,
        const T3* __restrict__ x, const ub32* __restrict__ n, const ub32* __restrict__ c, 
        const entryStructre_field<T,T3>* __restrict__ dataIn)
{
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_dir;
    ub32 th_id = threadIdx.x + (blockDim.x * blockIdx.x);
    T3 currentX = x[th_id];
    ub32 current_n;
            
    if(isPreRand)
    {
        current_n = n[th_id];
    }
    else
    {
        current_n = curand(state) % elemNum_smp_dir;
    }

    entryStructre_field<T,T3> pixelData = dataIn[current_n];

    thL[0].x = curr_constStr->aperture_kappa_l * pixelData.pixel.illuminationDir.x;
    thL[1].x = curr_constStr->aperture_kappa_l * pixelData.pixel.illuminationDir.y;
    thL[2].x = curr_constStr->aperture_kappa_l * pixelData.pixel.illuminationDir.z;

    thL[0].y = pixelData.pixel.k * (currentX.x - pixelData.pixel.illuminationP.x);
    thL[1].y = pixelData.pixel.k * (currentX.y - pixelData.pixel.illuminationP.y);
    thL[2].y = pixelData.pixel.k * (currentX.z - pixelData.pixel.illuminationP.z);

    dirv->x = pixelData.pixel.viewDir.x;
    dirv->y = pixelData.pixel.viewDir.y;
    dirv->z = pixelData.pixel.viewDir.z;

}

template<typename T, typename T3>
__global__ void sampleDirection_random(T3 *w0, T *pw0, ub64 seed)
{
    randomScatter scattering_structe = 0;
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    sample_direction<T,T3,3>(w0 + threadIdx.x, &scattering_structe, &state);
    pw0[threadIdx.x] = scattering_contribution<T,T3,3,false>(&scattering_structe, w0, w0);
}

template<typename T, typename T2, typename T3, typename correlationType>
__launch_bounds__(32)
__global__ void sampleDirection_gaussian_sum(T3* w0, T* pw0,
        const T3* __restrict__ x0, const ub32* __restrict__ n, const ub32* __restrict__ c,
        const correlationType* __restrict__ dataIn, ub64 seed)
{
    __shared__ T sampleProbabilitiesCdf[32][MAX_DIM*MAX_DIM];
    const constStructre<T>* curr_constStr = (const constStructre<T>*) &constStr_smp_dir;

    ub32 th_id = threadIdx.x + (blockDim.x * blockIdx.x);

    unsigned char mixtureIdx1, mixtureIdx2;
    T gamma_s;
    T2 gamma_s_over_beta_0;
    T2 thL[3];
    T3 real_conv_mu_1, real_conv_mu_2, real_conv_mu, w;
    T real_c_1, real_c_2, real_c, log_C;
    T2 C_tmp;
    T rkappa_1, rkappa_2, kappa;
    T alpha_sum = 0.0, alpha, curr_alpha_sum = 0.0;

    curandState_t state;
    curand_init(seed, th_id, 0, &state);

    T3 dirv;

    buildTh<T,T2,T3>(thL, &dirv, &state, x0, n, c, dataIn);

    __syncthreads();

    // calculate the convolution
    // no need to calculate attenuation, all mixtures have the same atternuation
    for(mixtureIdx1 = 0; mixtureIdx1 < curr_constStr->mixturesNum; mixtureIdx1++)
    {
        gamma_s = curr_constStr->mixtureMu[mixtureIdx1];
        if(abs(gamma_s) < 0.000000001)
        {
            gamma_s = 0.000000001;
        }
        
        gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
            complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
            complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
            complexSquare(cfma(gamma_s , dirv.z , thL[2])));

        real_conv_mu_1.x = realMult(gamma_s_over_beta_0,thL[0]);
        real_conv_mu_1.y = realMult(gamma_s_over_beta_0,thL[1]);
        real_conv_mu_1.z = realMult(gamma_s_over_beta_0,thL[2]);

        rkappa_1 = rnorm3d(real_conv_mu_1.x,real_conv_mu_1.y,real_conv_mu_1.z);

        w.x = real_conv_mu_1.x * rkappa_1;
        w.y = real_conv_mu_1.y * rkappa_1;
        w.z = real_conv_mu_1.z * rkappa_1;

        C_tmp = complexSqrt(
            complexSquare(cfma(gamma_s , w.x , thL[0])) +
            complexSquare(cfma(gamma_s , w.y , thL[1])) +
            complexSquare(cfma(gamma_s , w.z , thL[2])));

        real_c_1 = C_tmp.x - realComplexLog(C_tmp) + 
            curr_constStr->mixtureC[mixtureIdx1] - 
            fma(w.x , real_conv_mu_1.x , fma(w.y , real_conv_mu_1.y , w.z * real_conv_mu_1.z));

        for(mixtureIdx2 = 0; mixtureIdx2 < curr_constStr->mixturesNum; mixtureIdx2++)
        {
            gamma_s = curr_constStr->mixtureMu[mixtureIdx2];
            if(abs(gamma_s) < 0.000000001)
            {
                gamma_s = 0.000000001;
            }
        
            gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
                complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
                complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
                complexSquare(cfma(gamma_s , dirv.z , thL[2])));

            real_conv_mu_2.x = realMult(gamma_s_over_beta_0,thL[0]);
            real_conv_mu_2.y = realMult(gamma_s_over_beta_0,thL[1]);
            real_conv_mu_2.z = realMult(gamma_s_over_beta_0,thL[2]);

            rkappa_2 = rnorm3d(real_conv_mu_2.x,real_conv_mu_2.y,real_conv_mu_2.z);

            w.x = real_conv_mu_2.x * rkappa_2;
            w.y = real_conv_mu_2.y * rkappa_2;
            w.z = real_conv_mu_2.z * rkappa_2;

            C_tmp = complexSqrt(
                complexSquare(cfma(gamma_s , w.x , thL[0])) +
                complexSquare(cfma(gamma_s , w.y , thL[1])) +
                complexSquare(cfma(gamma_s , w.z , thL[2])));

            real_c_2 = C_tmp.x - realComplexLog(C_tmp) + 
                curr_constStr->mixtureC[mixtureIdx2] - 
                fma(w.x , real_conv_mu_2.x , fma(w.y , real_conv_mu_2.y , w.z * real_conv_mu_2.z));

            // multiple two mixtures
            real_conv_mu.x = real_conv_mu_1.x + real_conv_mu_2.x;
            real_conv_mu.y = real_conv_mu_1.y + real_conv_mu_2.y;
            real_conv_mu.z = real_conv_mu_1.z + real_conv_mu_2.z;
            real_c = real_c_1 + real_c_2;

            // cdf (not normalized) for each mixture
            kappa = norm3d(real_conv_mu.x,real_conv_mu.y,real_conv_mu.z);
            log_C = 0.5*log(kappa) - LOG_2_PIE_MULT_3_OVER_2 - logbesseli_05(kappa);

            alpha = exp(real_c - log_C - 2 * curr_constStr->aperture_kappa_l);
            alpha_sum += alpha;

            sampleProbabilitiesCdf[threadIdx.x][mixtureIdx2 + mixtureIdx1 * curr_constStr->mixturesNum] = alpha;
        }
    }

    // choose one of the samples
    T cdf_prob;
    randomUniform(&state,&cdf_prob);
    cdf_prob *= alpha_sum;
    
    unsigned short totalMixtures = curr_constStr->mixturesNum * curr_constStr->mixturesNum;

    for(unsigned short idxNum = 0; idxNum < totalMixtures ; idxNum++)
    {
        curr_alpha_sum += sampleProbabilitiesCdf[threadIdx.x][idxNum];
        if(curr_alpha_sum > cdf_prob)
        {
            mixtureIdx1 = idxNum / curr_constStr->mixturesNum;
            mixtureIdx2 = idxNum % curr_constStr->mixturesNum;
            break;
        }
    }

    // re-calculate normalized mu
            
    gamma_s = curr_constStr->mixtureMu[mixtureIdx1];
    if(abs(gamma_s) < 0.000000001)
    {
        gamma_s = 0.000000001;
	}
        
	gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
        complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
        complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
        complexSquare(cfma(gamma_s , dirv.z , thL[2])));

	real_conv_mu_1.x = realMult(gamma_s_over_beta_0,thL[0]);
    real_conv_mu_1.y = realMult(gamma_s_over_beta_0,thL[1]);
	real_conv_mu_1.z = realMult(gamma_s_over_beta_0,thL[2]);

    gamma_s = curr_constStr->mixtureMu[mixtureIdx2];
    if(abs(gamma_s) < 0.000000001)
    {
        gamma_s = 0.000000001;
	}
        
	gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
        complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
        complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
        complexSquare(cfma(gamma_s , dirv.z , thL[2])));

	real_conv_mu_2.x = realMult(gamma_s_over_beta_0,thL[0]);
    real_conv_mu_2.y = realMult(gamma_s_over_beta_0,thL[1]);
	real_conv_mu_2.z = realMult(gamma_s_over_beta_0,thL[2]);

    real_conv_mu.x = real_conv_mu_1.x + real_conv_mu_2.x;
    real_conv_mu.y = real_conv_mu_1.y + real_conv_mu_2.y;
    real_conv_mu.z = real_conv_mu_1.z + real_conv_mu_2.z;

    kappa = norm3d(real_conv_mu.x,real_conv_mu.y,real_conv_mu.z);

    real_conv_mu.x /= kappa;
    real_conv_mu.y /= kappa;
    real_conv_mu.z /= kappa;

    // sample a direction
    T3 sample_housed = random_vMF_direction<T,T2,T3>(real_conv_mu, kappa, &state);

    __syncthreads();

    // compute the probabilities
    T current_pw0 = 0.0;
    T muTimesX;

    for(mixtureIdx1 = 0; mixtureIdx1 < curr_constStr->mixturesNum; mixtureIdx1++)
    {
        gamma_s = curr_constStr->mixtureMu[mixtureIdx1];
        if(abs(gamma_s) < 0.000000001)
        {
            gamma_s = 0.000000001;
        }
        
        gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
            complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
            complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
            complexSquare(cfma(gamma_s , dirv.z , thL[2])));

        real_conv_mu_1.x = realMult(gamma_s_over_beta_0,thL[0]);
        real_conv_mu_1.y = realMult(gamma_s_over_beta_0,thL[1]);
        real_conv_mu_1.z = realMult(gamma_s_over_beta_0,thL[2]);

        for(mixtureIdx2 = 0; mixtureIdx2 < curr_constStr->mixturesNum; mixtureIdx2++)
        {
            gamma_s = curr_constStr->mixtureMu[mixtureIdx2];
            if(abs(gamma_s) < 0.000000001)
            {
                gamma_s = 0.000000001;
            }
        
            gamma_s_over_beta_0 = gamma_s * rComplexSqrt(
                complexSquare(cfma(gamma_s , dirv.x , thL[0])) +
                complexSquare(cfma(gamma_s , dirv.y , thL[1])) +
                complexSquare(cfma(gamma_s , dirv.z , thL[2])));

            real_conv_mu_2.x = realMult(gamma_s_over_beta_0,thL[0]);
            real_conv_mu_2.y = realMult(gamma_s_over_beta_0,thL[1]);
            real_conv_mu_2.z = realMult(gamma_s_over_beta_0,thL[2]);

            // multiple two mixtures
            real_conv_mu.x = real_conv_mu_1.x + real_conv_mu_2.x;
            real_conv_mu.y = real_conv_mu_1.y + real_conv_mu_2.y;
            real_conv_mu.z = real_conv_mu_1.z + real_conv_mu_2.z;

            // cdf (not normalized) for each mixture
            kappa = norm3d(real_conv_mu.x,real_conv_mu.y,real_conv_mu.z);
            log_C = 0.5*log(kappa) - LOG_2_PIE_MULT_3_OVER_2 - logbesseli_05(kappa);
            muTimesX = real_conv_mu.x * sample_housed.x + 
                       real_conv_mu.y * sample_housed.y +
                       real_conv_mu.z * sample_housed.z;

            current_pw0 += (sampleProbabilitiesCdf[threadIdx.x][mixtureIdx2 + mixtureIdx1 * curr_constStr->mixturesNum] / alpha_sum) *
                exp(muTimesX + log_C);
        }
    }

    w0[th_id] = sample_housed;
    pw0[th_id] = current_pw0;
}

template<typename T, typename T2, typename T3>
void sampleDirection(const T3 *x0, T3 *w0, T* pw0, const ub32 *n, const ub32 *c,
        ub32 is_correlation, ub32 sample_direction_flag,
        const void* globalMem, ub64 curr_seed, T min_prob)
{
    if(sample_direction_flag == 1)
    {
        sampleDirection_random<T,T3><<<1, THREADS_NUM>>>(w0,pw0,curr_seed);
    }
    else if(sample_direction_flag == 4)
    {
        if(is_correlation)
        {
            entryStructre_correlation<T,T3> *globalMem_corr = (entryStructre_correlation<T,T3>*) globalMem;

            sampleDirection_gaussian_sum<T,T2,T3,entryStructre_correlation<T,T3>><<<32, THREADS_NUM/32>>>(
                    w0,pw0,x0,n,c,globalMem_corr,curr_seed);
        }
        else
        {
            entryStructre_field<T,T3> *globalMem_field = (entryStructre_field<T,T3>*) globalMem;

            sampleDirection_gaussian_sum<T,T2,T3,entryStructre_field<T,T3>><<<32, THREADS_NUM/32>>>(
                    w0,pw0,x0,n,c,globalMem_field,curr_seed);
        }
    }

    if(!is_correlation)
    {
        // in field we take the sqrt of pw0
        global_sqrt<T><<<1, THREADS_NUM>>>(pw0);
    }

    global_min<T><<<1, THREADS_NUM>>>(pw0, min_prob);
}

template<typename T>
void sampleDirection_setConstMem(constStructre<T> *constMem, ub32 number_of_elements, bool ipr)
{
    cudaMemcpyToSymbol(constStr_smp_dir, constMem, sizeof(constStructre<T>));
    cudaMemcpyToSymbol(elemNum_smp_dir, &number_of_elements, sizeof(ub32));
    cudaMemcpyToSymbol(isPreRand, &ipr, sizeof(bool));
}

template void sampleDirection<double,double2,double3>(const double3 *x0, double3 *w0, double* pw0, const ub32 *n, const ub32 *c,
        ub32 is_correlation, ub32 sample_direction_flag, const void* globalMem, ub64 curr_seed, double min_prob);

template void sampleDirection<float,float2,float3>(const float3 *x0, float3 *w0, float* pw0, const ub32 *n, const ub32 *c,
        ub32 is_correlation, ub32 sample_direction_flag, const void* globalMem, ub64 curr_seed, float min_prob);

template void sampleDirection_setConstMem<double>(constStructre<double> *constMem, ub32 number_of_elements, bool ipr);
template void sampleDirection_setConstMem<float>(constStructre<float> *constMem, ub32 number_of_elements, bool ipr);
