#include "tabulated_scattering.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

__constant__ ub32 tab_entries_cons;
__constant__ ub32 tab_rounds_cons; // ceil(tab_entries_cons / 1024)
 
__device__ ub32 index_in_between(int n)
{
    if(n < 0)
    {
        n = 0;
    }

    if(n >= tab_entries_cons)
    {
        n = tab_entries_cons - 1;
    }

    return (ub32) n;
}
        
__device__ ub32 round_num(double a)
{
    return index_in_between(__double2int_rn(a));
}

__device__ ub32 round_num(float a)
{
    return index_in_between(__float2int_rn(a));
}
        
template<typename T, typename T2, ub32 dims>
__global__ void init_tabulated_pdf_cdf(tabulatedScatter<T>* tabScatterer, T* tmp_cdf, const T2* tabulated_values)
{
    __shared__ volatile T pdf_vals[THREADS_NUM];
    __shared__ volatile T cdf_vals[THREADS_NUM];
    __shared__ T cdf_sum_vals[THREADS_NUM];
    __shared__ T cdf_total_sum;
    __shared__ volatile T cdf_tmp_sum[32];

    for(ub32 blockNum = 0; blockNum < tab_rounds_cons; blockNum++)
    {
        ub32 transpose_thread = (threadIdx.x % 32) * 32 + threadIdx.x / 32;
        ub32 currentSample = threadIdx.x + blockNum * THREADS_NUM;
        if(currentSample < tab_entries_cons)
        {
            T2 curr_value = tabulated_values[currentSample];
            if(dims == 2)
            {
                pdf_vals[transpose_thread] = curr_value.x * curr_value.x + curr_value.y * curr_value.y;
            }
            if(dims == 3)
            {
                T sintheta = sinpi(((T)currentSample) / (tab_entries_cons - 1));
                pdf_vals[transpose_thread] = (curr_value.x * curr_value.x + curr_value.y * curr_value.y) * sintheta;
            }
        }
        else
        {
            pdf_vals[transpose_thread] = 0;
        }
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
                if(blockNum == 0) {cdf_total_sum = (T) 0.0;}
                cdf_tmp_sum[0] = (T) 0.0;
                for(ub32 tmpCdfSumIter = 1; tmpCdfSumIter < 32; tmpCdfSumIter++ )
                {
                    cdf_tmp_sum[tmpCdfSumIter] = cdf_tmp_sum[tmpCdfSumIter - 1] + cdf_vals[991 + tmpCdfSumIter];
                }
                cdf_sum_vals[blockNum] = cdf_vals[1023] + cdf_tmp_sum[31];
                cdf_total_sum += cdf_sum_vals[blockNum];
            }    
        }
        __syncthreads();

        if(currentSample < tab_entries_cons)
        {
            T2 curr_value = tabulated_values[currentSample];
            tmp_cdf[currentSample] = cdf_vals[transpose_thread] + cdf_tmp_sum[threadIdx.x / 32];
            tabScatterer->evalAmp[currentSample] = curr_value.x * curr_value.x + curr_value.y * curr_value.y;
        }
    }

    __syncthreads();
    T addCdfVal = (T)0.0;
    T cs = (cdf_total_sum / ((T) tab_entries_cons)) * 2 * M_PI;

    if(dims == 3)
    {
        cs *= M_PI;
    }
            
    for(ub32 blockNum = 0; blockNum < tab_rounds_cons; blockNum++)
    {
        ub32 currentSample = threadIdx.x + blockNum * THREADS_NUM;
        if(currentSample < tab_entries_cons)
        {
            tabScatterer->evalAmp[currentSample] /= cs;
            tmp_cdf[currentSample] = (addCdfVal + tmp_cdf[currentSample] ) / cdf_total_sum;
        }
        addCdfVal += cdf_sum_vals[blockNum];
    }
    
}

template<typename T>
__global__ void init_tabulated_icdf(tabulatedScatter<T>* tabScatterer, T* tmp_cdf, bool is3D)
{
    ub32 currentSample = threadIdx.x + blockIdx.x * THREADS_NUM;

    if(currentSample < tab_entries_cons)
    {
        T x_axis = currentSample / ((T)(tab_entries_cons - 1)); // from 0 to 1 (including)
                
        // for y_axis, we search x_axis value in tmp_cdf index, using binary search
        ub32 y_axis_idx = binary_search<T>(tmp_cdf, x_axis, tab_entries_cons);

        // and the value of y_axis is from [0,pi] or [0,2*pi]
        T y_axis = (((T)(y_axis_idx)) / (tab_entries_cons - 1)) * (is3D ? M_PI : (2*M_PI));

        // compute the cos of iCdf, because that what we need for sampling
        tabScatterer->sampleIcdf[currentSample] = cos(y_axis);
    }
}

template<typename T, typename T2>
void* init_tabulated_scatterer(ub32 tabulated_entries_num, const T2* tabulated_values, bool is3D)
{
    tabulatedScatter<T> tabScatterer;
    tabulatedScatter<T> *gpu_tabScatterer;
    T* tmp_cdf;
    ub32 tab_rounds = (tabulated_entries_num - 1) /THREADS_NUM + 1;
    T2 *tabulated_values_gpu;
    cudaMalloc(&tabulated_values_gpu, tabulated_entries_num * sizeof(T2));
    cudaMemcpy(tabulated_values_gpu, tabulated_values, tabulated_entries_num * sizeof(T2), cudaMemcpyHostToDevice);

    tabScatterer.tabulated_entries = tabulated_entries_num;
    cudaMalloc(&tabScatterer.evalAmp, tabulated_entries_num * sizeof(T));
    cudaMalloc(&tabScatterer.sampleIcdf, tabulated_entries_num * sizeof(T));
    cudaMalloc(&tmp_cdf, tabulated_entries_num * sizeof(T));

    cudaMalloc(&gpu_tabScatterer, sizeof(tabulatedScatter<T>));
    cudaMemcpy(gpu_tabScatterer, &tabScatterer, sizeof(tabulatedScatter<T>), cudaMemcpyHostToDevice);

    cudaMemcpyToSymbol(tab_entries_cons, &tabulated_entries_num, sizeof(ub32));
    cudaMemcpyToSymbol(tab_rounds_cons, &tab_rounds, sizeof(ub32));

    if(is3D)
    {
        init_tabulated_pdf_cdf<T,T2,3><<<1,THREADS_NUM>>>(gpu_tabScatterer, tmp_cdf, tabulated_values_gpu);
    }
    else
    {
        init_tabulated_pdf_cdf<T,T2,2><<<1,THREADS_NUM>>>(gpu_tabScatterer, tmp_cdf, tabulated_values_gpu);
    }

    init_tabulated_icdf<T><<<tab_rounds,THREADS_NUM>>>(gpu_tabScatterer, tmp_cdf, is3D);
    
    cudaFree(tabulated_values_gpu);
    cudaFree(tmp_cdf);
    return (void *) gpu_tabScatterer;
}

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const tabulatedScatter<T>* scattering_struct, const T3* D, const T3* W)
{
    T ret_val;
    if(dims == 3)
    {
        T cosang = fma(D->x , W->x ,
                   fma(D->y , W->y ,
                       D->z * W->z));

        T ang = acos(cosang); // scat_angle
        ret_val = scattering_struct->evalAmp[round_num(ang * (tab_entries_cons - 1) / (M_PI))];
    }
    if(dims == 2)
    {
        T cosang = fma(D->x , W->x ,
                      (D->y * W->y));

        T ang = acos(cosang); // scat_angle
        ret_val = scattering_struct->evalAmp[round_num(ang * (tab_entries_cons - 1) / (2*M_PI))];
    }

    if(is_sqrt)
    {
        ret_val = sqrt(ret_val);
    }

    return ret_val;
}

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const tabulatedScatter<T>* sample_function, ub64 seed, ub32 dims_num)
{
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    if(dims_num == 3)
    {
        sample_direction<T,T3,3>(w + threadIdx.x, sample_function, &state);
    }
    if(dims_num == 2)
    {
        sample_direction<T,T3,2>(w + threadIdx.x, sample_function, &state);
    }
    
}

template<typename T, typename T3, ub32 dims>
__device__ void sample_direction(T3 *w, const tabulatedScatter<T>* sample_function, curandState_t *state)
{
    T random_num;
    randomUniform(state, &random_num);

    T costheta = sample_function->sampleIcdf[round_num(random_num * (tab_entries_cons - 1))];

    rotate_by_theta<T,T3,dims>(w, costheta, state);
}

template void* init_tabulated_scatterer<double,double2>(ub32 tabulated_entries_num, const double2* tabulated_values, bool is3D);
template void* init_tabulated_scatterer<float,float2>(ub32 tabulated_entries_num, const float2* tabulated_values, bool is3D);

template __device__ double scattering_contribution<double,double3,2,true>(const tabulatedScatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,true>(const tabulatedScatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,true>(const tabulatedScatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,true>(const tabulatedScatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,2,false>(const tabulatedScatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,false>(const tabulatedScatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,false>(const tabulatedScatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,false>(const tabulatedScatter<float>* scattering_struct, const float3* D, const float3* W);

template __global__ void sample_direction<double,double3>(double3 *w, const tabulatedScatter<double>* sample_function, ub64 seed, ub32 dims);
template __global__ void sample_direction<float,float3>(float3 *w, const tabulatedScatter<float>* sample_function, ub64 seed, ub32 dims);

template __device__ void sample_direction<double,double3,2>(double3 *w, const tabulatedScatter<double>* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,2>(float3 *w, const tabulatedScatter<float>* sample_function, curandState_t *state);
template __device__ void sample_direction<double,double3,3>(double3 *w, const tabulatedScatter<double>* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,3>(float3 *w, const tabulatedScatter<float>* sample_function, curandState_t *state);
