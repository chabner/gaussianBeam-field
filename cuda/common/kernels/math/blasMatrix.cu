#include "blasMatrix.h"

// compute C = A * B;
// where A in size of [n x m] and B in size of [m x k]
// A and B are stored as vectors, in such way that:
// A[0]:i_n = 0, i_m = 0. A[1]:i_n = 0, i_m = 1, ... A[m]:i_n = 1, i_m = 0

void matrixMult(cublasHandle_t cublas_handle, ub32 n, ub32 m, ub32 k,
	const float2* A, const float2* B, float2* C)
{
    cuComplex alpha = make_cuComplex(1.0f,0.0f);
    cuComplex beta  = make_cuComplex(0.0f,0.0f);

    cublasCgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
        k, n, m, &alpha, (cuComplex*) B, k, (cuComplex*) A, m, &beta,
        (cuComplex*) C, k); 
}

void matrixMult(cublasHandle_t cublas_handle, ub32 n, ub32 m, ub32 k,
	const double2* A, const double2* B, double2* C)
{
    cuDoubleComplex alpha = make_cuDoubleComplex(1.0,0.0);
    cuDoubleComplex beta  = make_cuDoubleComplex(0.0,0.0);

    cublasZgemm(cublas_handle, CUBLAS_OP_N, CUBLAS_OP_N,
        k, n, m, &alpha, (cuDoubleComplex*) B, k, (cuDoubleComplex*) A, m, &beta,
        (cuDoubleComplex*) C, k); 
}

// compute A = A + B;
// where A and B are vectors in size of n
void matrixAdd(cublasHandle_t cublas_handle, ub32 n,
	const float2* alpha, float2* A, const float2* B)
{
    cublasCaxpy(cublas_handle, n, (cuComplex*) alpha,
        (cuComplex*) B, 1,(cuComplex*) A, 1);
}

void matrixAdd(cublasHandle_t cublas_handle, ub32 n,
	const double2* alpha, double2* A, const double2* B)
{
    cublasZaxpy(cublas_handle, n, (cuDoubleComplex*) alpha,
        (cuDoubleComplex*) B, 1,(cuDoubleComplex*) A, 1);
}
