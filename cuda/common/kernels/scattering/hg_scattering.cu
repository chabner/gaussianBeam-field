#include "hg_scattering.h"
#include "cuMath.cu"
#include "cuRandoms.cu"

template<typename T, typename T3, ub32 dims, bool is_sqrt>
__device__ T scattering_contribution(const HGscatter<T>* scattering_struct, const T3* D, const T3* W)
{
    T ret_val;
    if(dims == 3)
    {
        T cosang = fma(D->x , W->x ,
                   fma(D->y , W->y ,
                       D->z * W->z));

        T downCalc = fma(scattering_struct->g_down2, cosang, scattering_struct->g_down1);

        ret_val = scattering_struct->g_up / (downCalc * sqrt(downCalc));
    }
    if(dims == 2)
    {
        T cosang = fma(D->x , W->x ,
                      (D->y * W->y));

        T downCalc = fma(scattering_struct->g_down2, cosang, scattering_struct->g_down1);

        ret_val = (2.0 * scattering_struct->g_up) / (downCalc);
    }

    if(is_sqrt)
    {
        ret_val = sqrt(ret_val);
    }

    return ret_val;
}

template<typename T>
void* init_hg_scatterer(T g)
{
    HGscatter<T> hgScatterer;
    void *gpu_hgScatterer = 0;

    hgScatterer.g = g;
    hgScatterer.g_up = (T) ((1.0 - g * g) / (4.0 * M_PI));
    hgScatterer.g_down1 = (T) (1.0 + g * g);
    hgScatterer.g_down2 = (T) (-2.0 * g);

    cudaMalloc(&gpu_hgScatterer, sizeof(HGscatter<T>));
    cudaMemcpy(gpu_hgScatterer, &hgScatterer, sizeof(HGscatter<T>), cudaMemcpyHostToDevice);

    return gpu_hgScatterer;
}

template<typename T, typename T3>
__global__ void sample_direction(T3 *w, const HGscatter<T>* sample_function, ub64 seed, ub32 dims_num)
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
__device__ void sample_direction(T3 *w, const HGscatter<T>* sample_function, curandState_t *state)
{
    T g = sample_function->g;

    if(abs(g) < 0.001)
    {
        if(dims == 3)
        {
            *w = randomDirection3<T,T3>(state);
        }

        if(dims == 2)
        {
            *w = randomDirection2<T,T3>(state);
        }
    }
    else
    {
        T random_num, costheta;
        randomUniform(state, &random_num);
        if(dims == 3)
        {
            T s = (1.0 - g * g) / (1.0 - g + 2 * g * random_num);
            costheta = (1.0 + g * g - s * s) / (2.0 * g);
        }
        if(dims == 2)
        {
            costheta = cos(2.0*atan( ((1.0-g)/(1.0+g)) * tan( (M_PI/2.0) * (1-2.0*random_num)) ));
        }

        rotate_by_theta<T,T3,dims>(w, costheta, state);
    }
}

template __device__ double scattering_contribution<double,double3,2,false>(const HGscatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,false>(const HGscatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,false>(const HGscatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,false>(const HGscatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,2,true>(const HGscatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,2,true>(const HGscatter<float>* scattering_struct, const float3* D, const float3* W);
template __device__ double scattering_contribution<double,double3,3,true>(const HGscatter<double>* scattering_struct, const double3* D, const double3* W);
template __device__ float scattering_contribution<float,float3,3,true>(const HGscatter<float>* scattering_struct, const float3* D, const float3* W);

template void* init_hg_scatterer<double>(double g);
template void* init_hg_scatterer<float>(float g);

template __global__ void sample_direction<double,double3>(double3 *w, const HGscatter<double>* sample_function, ub64 seed, ub32 dims);
template __global__ void sample_direction<float,float3>(float3 *w, const HGscatter<float>* sample_function, ub64 seed, ub32 dims);

template __device__ void sample_direction<double,double3,2>(double3 *w, const HGscatter<double>* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,2>(float3 *w, const HGscatter<float>* sample_function, curandState_t *state);
template __device__ void sample_direction<double,double3,3>(double3 *w, const HGscatter<double>* sample_function, curandState_t *state);
template __device__ void sample_direction<float,float3,3>(float3 *w, const HGscatter<float>* sample_function, curandState_t *state);
