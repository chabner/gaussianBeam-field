#pragma once

#include "standardTypes.h"
#include "cublas_v2.h"

// compute C = A * B;
// where A in size of [n x m] and B in size of [m x k]
// A and B are stored as vectors, in such way that:
// A[0]:i_n = 0, i_m = 0. A[1]:i_n = 0, i_m = 1, ... A[m]:i_n = 1, i_m = 0
void matrixMult(cublasHandle_t cublas_handle, ub32 n, ub32 m, ub32 k, 
        const float2* A, const float2* B, float2* C);
void matrixMult(cublasHandle_t cublas_handle, ub32 n, ub32 m, ub32 k, 
        const double2* A, const double2* B, double2* C);

// compute A = A + alpha * B;
// where A and B are vectors in size of n
void matrixAdd(cublasHandle_t cublas_handle, ub32 n,
        const float2* alpha, float2* A, const float2* B);
void matrixAdd(cublasHandle_t cublas_handle, ub32 n,
        const double2* alpha, double2* A, const double2* B);
