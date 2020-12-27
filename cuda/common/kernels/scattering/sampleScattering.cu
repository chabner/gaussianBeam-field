#include "sampleScattering.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

// Random Phase

template<typename T, typename T2>
__global__ void const_path_contribution_corr(T2 *const_path, const T* px0)
{
    const_path[threadIdx.x] = buildComplex(1/px0[threadIdx.x]);
}

template<typename T, typename T2>
__global__ void const_path_contribution_field(T2 *const_path, const T* px0, ub64 seed)
{
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    T randNum;
    T zero = 0.0;
    randomUniform(&state, &randNum);
    randNum = randNum * 2 * M_PI;

    const_path[threadIdx.x] = (1/px0[threadIdx.x]) * complexExponent(buildComplex(zero,randNum));
}

template<typename T, typename T2>
void const_path_contribution(T2 *const_path, const T* px0, ub32 is_correlation, ub64 seed)
{
    if(is_correlation)
    {
        const_path_contribution_corr<T,T2><<<1, THREADS_NUM>>>(const_path, px0);
    }
    else
    {
        const_path_contribution_field<T,T2><<<1, THREADS_NUM>>>(const_path, px0, seed);
    }
}

template<typename T>
void px0_mult_pw0(T *px0, const T* pw0)
{
    global_mult<T><<<1, THREADS_NUM>>>(px0, pw0);
}

// x scattering propagate
        
template<typename T, typename T3, ub32 dims_num>
__global__ void propagate_x_kernel(T3 *x, pathsList<T> *pl, const T3* w, T sigt, ub64 seed)
{
    curandState_t state;
    curand_init(seed, threadIdx.x, 0, &state);

    T randNum;
    randomUniform(&state, &randNum);

    randNum = fmin(randNum,(T)0.999);

    T d = -log(-randNum + 1.0)/(sigt);

    x[threadIdx.x].x += d * w[threadIdx.x].x;
    x[threadIdx.x].y += d * w[threadIdx.x].y;
    if(dims_num == 3)
    {
        x[threadIdx.x].z += d * w[threadIdx.x].z;
    }

    pl->path_length[threadIdx.x] += d;
}
        
template<typename T, typename T3>
void propagate_x(T3 *x, pathsList<T> *pl, const T3* w, T sigt, ub64 seed, ub32 dims_num)
{
    if(dims_num == 2)
    {
        propagate_x_kernel<T,T3,2><<<1, THREADS_NUM>>>(x, pl, w, sigt, seed);
    }

    if(dims_num == 3)
    {
        propagate_x_kernel<T,T3,3><<<1, THREADS_NUM>>>(x, pl, w, sigt, seed);
    }
}
 
// activated paths
template<typename T>
__global__ void initiate_path_list_kernel(pathsList<T> *pl)
{
    pl->paths_inside[threadIdx.x] = 0xFFFFFFFF;
}

template<typename T>
void initiate_path_list(pathsList<T> *pl)
{
    initiate_path_list_kernel<T><<<1, 32>>>(pl);
    cudaMemset(pl->path_length, (T) 0.0, THREADS_NUM * sizeof(T));
}

template<typename T, typename T3, ub32 dims_num>
__global__ void is_inside_path_list_kernel(pathsList<T> *pl, const T3 *x, T3 box_min, T3 box_max)
{
    __shared__ volatile ub32 paths_sum[32];
    __shared__ ub32 is_inside[THREADS_NUM];

    ub32 bit_num = (threadIdx.x % 32);
    ub32 p_num = threadIdx.x / 32;

    T3 curr_x = x[threadIdx.x];

    if(pl->paths_inside[p_num])
    {
        bool check_if_inside;
        if(dims_num == 3)
        {
            check_if_inside = (curr_x.x >= box_min.x) && (curr_x.x <= box_max.x) &&
                              (curr_x.y >= box_min.y) && (curr_x.y <= box_max.y) &&
                              (curr_x.z >= box_min.z) && (curr_x.z <= box_max.z);
        }
        if(dims_num == 2)
        {
            check_if_inside = (curr_x.x >= box_min.x) && (curr_x.x <= box_max.x) &&
                              (curr_x.y >= box_min.y) && (curr_x.y <= box_max.y);
        }

        is_inside[threadIdx.x] = check_if_inside << bit_num;
    }
    else
    {
        is_inside[threadIdx.x] = 0;
    }
    __syncthreads();

    if(threadIdx.x < 32)
    {
        ub32 thread_paths_inside = 0;
        ub32 thread_sq = 32 * threadIdx.x;

        thread_paths_inside += is_inside[     thread_sq];
        thread_paths_inside += is_inside[1  + thread_sq];
        thread_paths_inside += is_inside[2  + thread_sq];
        thread_paths_inside += is_inside[3  + thread_sq];
        thread_paths_inside += is_inside[4  + thread_sq];
        thread_paths_inside += is_inside[5  + thread_sq];
        thread_paths_inside += is_inside[6  + thread_sq];
        thread_paths_inside += is_inside[7  + thread_sq];
        thread_paths_inside += is_inside[8  + thread_sq];
        thread_paths_inside += is_inside[9  + thread_sq];
        thread_paths_inside += is_inside[10 + thread_sq];
        thread_paths_inside += is_inside[11 + thread_sq];
        thread_paths_inside += is_inside[12 + thread_sq];
        thread_paths_inside += is_inside[13 + thread_sq];
        thread_paths_inside += is_inside[14 + thread_sq];
        thread_paths_inside += is_inside[15 + thread_sq];
        thread_paths_inside += is_inside[16 + thread_sq];
        thread_paths_inside += is_inside[17 + thread_sq];
        thread_paths_inside += is_inside[18 + thread_sq];
        thread_paths_inside += is_inside[19 + thread_sq];
        thread_paths_inside += is_inside[20 + thread_sq];
        thread_paths_inside += is_inside[21 + thread_sq];
        thread_paths_inside += is_inside[22 + thread_sq];
        thread_paths_inside += is_inside[23 + thread_sq];
        thread_paths_inside += is_inside[24 + thread_sq];
        thread_paths_inside += is_inside[25 + thread_sq];
        thread_paths_inside += is_inside[26 + thread_sq];
        thread_paths_inside += is_inside[27 + thread_sq];
        thread_paths_inside += is_inside[28 + thread_sq];
        thread_paths_inside += is_inside[29 + thread_sq];
        thread_paths_inside += is_inside[30 + thread_sq];
        thread_paths_inside += is_inside[31 + thread_sq];

        thread_paths_inside &= pl->paths_inside[threadIdx.x];
        pl->paths_inside[threadIdx.x] = thread_paths_inside;

        paths_sum[threadIdx.x] = __popc(thread_paths_inside);

        if(threadIdx.x >= 1 ) {paths_sum[threadIdx.x] += paths_sum[threadIdx.x - 1 ]; }
        if(threadIdx.x >= 2 ) {paths_sum[threadIdx.x] += paths_sum[threadIdx.x - 2 ]; }
        if(threadIdx.x >= 4 ) {paths_sum[threadIdx.x] += paths_sum[threadIdx.x - 4 ]; }
        if(threadIdx.x >= 8 ) {paths_sum[threadIdx.x] += paths_sum[threadIdx.x - 8 ]; }
        if(threadIdx.x >= 16) {paths_sum[threadIdx.x] += paths_sum[threadIdx.x - 16]; }
        if(threadIdx.x == 31) {pl->path_count = paths_sum[31];}
    }
    __syncthreads();

    ub32 current_paths_sum = (p_num == 0 ? 0 : paths_sum[p_num - 1]);
    ub32 thread_bit_mask = (1U << bit_num);

    ub32 thread_inside_paths = pl->paths_inside[p_num];

    // if the path which numbered as threadIdx.x is still active
    if (thread_inside_paths & thread_bit_mask)
    {
        // then put it in the list of active threads
        pl->paths_numbers_list[current_paths_sum + __popc((thread_bit_mask - 1) & thread_inside_paths)] = threadIdx.x;
    }
}

template<typename T, typename T3>
ub32 is_inside_path_list(pathsList<T> *pl, const T3 *x, T3 box_min, T3 box_max, ub32 dims_num)
{
    ub32 path_count;
    if(dims_num == 2)
    {
        is_inside_path_list_kernel<T,T3,2><<<1, THREADS_NUM>>>(pl, x, box_min, box_max);
    }
    if(dims_num == 3)
    {
        is_inside_path_list_kernel<T,T3,3><<<1, THREADS_NUM>>>(pl, x, box_min, box_max);
    }
    cudaMemcpy(&path_count, &pl->path_count, sizeof(ub32), cudaMemcpyDeviceToHost);

    return path_count;
}

// scattering function
template<typename T, typename T3>
void sample_new_directions(T3 *w, ub32 sample_flag, const void *sample_function, ub64 curr_seed, ub32 dims_num )
{
    if(sample_flag == 1)
    {
        sample_direction<T,T3><<<1, THREADS_NUM>>>(w, (const randomScatter*) sample_function, curr_seed, dims_num);
    }
    else if(sample_flag == 2)
    {
        sample_direction<T,T3><<<1, THREADS_NUM>>>(w, (const tabulatedScatter<T>*) sample_function, curr_seed, dims_num );
    }
    else if(sample_flag == 3)
    {
        sample_direction<T,T3><<<1, THREADS_NUM>>>(w, (const HGscatter<T>*) sample_function, curr_seed, dims_num );
    }


    if(dims_num == 2)
    {
        normalize_dir<T,T3,2><<<1, THREADS_NUM>>>(w);
    }
    if(dims_num == 3)
    {
        normalize_dir<T,T3,3><<<1, THREADS_NUM>>>(w);
    }
}

// explicit declarations
template void const_path_contribution<double,double2>(double2 *const_path, const double* px0, ub32 is_correlation, ub64 curr_seed);
template void const_path_contribution<float,float2>(float2 *const_path, const float* px0, ub32 is_correlation, ub64 curr_seed);

template void propagate_x<double,double3>(double3 *x, pathsList<double> *pl ,const double3* w, double sigt, ub64 seed, ub32 dims_num);
template void propagate_x<float,float3>(float3 *x, pathsList<float> *pl, const float3* w, float sigt, ub64 seed, ub32 dims_num);

template ub32 is_inside_path_list<double,double3>(pathsList<double> *pl, const double3 *x, double3 box_min, double3 box_max, ub32 dims_num);
template ub32 is_inside_path_list<float,float3>(pathsList<float> *pl, const float3 *x, float3 box_min, float3 box_max, ub32 dims_num);

template void sample_new_directions<double,double3>(double3 *w, ub32 sample_flag, const void *sample_function, ub64 curr_seed, ub32 dims_num );
template void sample_new_directions<float,float3>(float3 *w, ub32 sample_flag, const void *sample_function, ub64 curr_seed, ub32 dims_num );

template void px0_mult_pw0<double>(double *px0, const double* pw0);
template void px0_mult_pw0<float>(float *px0, const float* pw0);

template void initiate_path_list<double>(pathsList<double> *pl);
template void initiate_path_list<float>(pathsList<float> *pl);
