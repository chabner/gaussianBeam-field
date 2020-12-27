#pragma once

#include<stdint.h>

#include "cuda_runtime.h"

typedef uint16_t ub16;
typedef uint32_t ub32;
typedef uint64_t ub64;

#define THREADS_NUM 1024

#define LOG_2_PI 1.837877066409345483
#define M_PI 3.14159265358979323846
#define SQRT_2_OVER_PI 0.79788456080286535588
#define SQRT_2_OVER_PI_OVER2 0.39894228040143267794
#define LOG_2_PIE_MULT_3_OVER_2 2.75681559961401822534
