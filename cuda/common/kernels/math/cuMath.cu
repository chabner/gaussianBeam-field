#pragma once

#include "math.h"
#include "cuda_runtime.h"
#include "standardTypes.h"

// DOUBLES

__inline __device__ double2 buildComplex(double a, double b)
{
    double2 res;
    res.x = a; res.y = b;
    return res;
}

__inline __device__ double2 buildComplex(double a)
{
    return buildComplex(a,0.0);
}

__inline __device__ double2 imagAdd(double2 a, double b)
{
    a.y += b;
    return a;
}

__inline __device__ double2 operator +(double2 a, double2 b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

__inline __device__ double2 operator -(double2 a, double2 b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

__inline __device__ double2 operator -(double a, double2 b)
{
    b.x = a - b.x;
    b.y = -b.y;
    return b;
}

__inline __device__ double2 operator +(double2 a, double b)
{
    a.x += b;
    return a;
}

__inline__ __device__ double2 operator /(double2 a, double2 b)
{
    double2 c;
    double denominator = 1 / (fma(b.x, b.x, b.y * b.y));
    c.x = (fma(a.x, b.x, a.y * b.y)) * denominator;
    c.y = (fma(a.y, b.x, -a.x * b.y)) * denominator;
        
//     c.x = (fma(a.x, b.x, a.y * b.y)) / (fma(b.x, b.x, b.y * b.y));
//     c.y = (fma(a.y, b.x, -a.x * b.y)) / (fma(b.x, b.x, b.y * b.y));

    return c;
}

__inline __device__ double2 operator /(double2 a, double b)
{
    a.x /= b;
    a.y /= b;
    return a;
}


__inline__ __device__ double2 operator /(double a, double2 b)
{
    a /= (fma(b.x, b.x, b.y * b.y));
    b.x *= a;
    b.y *= (-a);
    return b;
}


__inline__ __device__ double2 complexr(double2 b)
{
    double a = 1.0/(fma(b.x, b.x, b.y * b.y));
    b.x *= a;
    b.y *= (-a);
    return b;
}


__inline __device__ double2 operator *(double a, double2 b)
{
    b.x *= a;
    b.y *= a;
    return b;
}


__inline __device__ double2 operator *(double2 a, double2 b)
{
    double2 tmpVar = a;

    a.x = fma(a.x, b.x, -a.y * b.y);
    a.y = fma(tmpVar.x, b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ double realMult(double2 a, double2 b)
{
    return fma(a.x, b.x, -a.y * b.y);
}

__inline__  __device__ double2 conjMult(double2 a, double2 b)
{
    double2 tmpVar = a;

    a.x = fma(a.x, b.x, a.y * b.y);
    a.y = fma(tmpVar.x, -b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ double2 complexConj(double2 a)
{
    a.y *= -1;
    return a;
}

__inline __device__ double2 complexSquare(double2 a)
{
    double2 res;
    res.x = (a.x - a.y) * (a.x + a.y);
//     res.x = fma(a.x,a.x,-a.y*a.y);
    res.y = 2 * a.x * a.y;
    return res;
}

__inline__ __device__ double2 complexSqrt(double2 a)
{
    double r, absrz;
    double2 res;

    if(a.y == 0)
    {
        res.x = (a.x >= 0 ? sqrt(a.x) : 0);
        res.y = (a.x >= 0 ? 0 : sqrt(-a.x));
    }
    else
    {
        r = sqrt(fma(a.x, a.x, a.y * a.y));
        a.x = a.x + r;
    //     absrz = 1 / sqrt((fma(r + a.x, r + a.x, a.y * a.y)) / r);
        absrz = 1 / sqrt((fma(a.x, a.x, a.y * a.y)) / r);
    //     res.x = absrz * (a.x + r);
        res.x = absrz * a.x;
        res.y = absrz * a.y;
    }
    return res;
}

__inline__ __device__ double2 rComplexSqrt(double2 a)
{
    double r, absrz;
    double2 res;
        
    if(a.y == 0)
    {
        res.x = (a.x >= 0 ? 1/sqrt(a.x) : 0);
        res.y = (a.x >= 0 ? 0 : 1/sqrt(-a.x));
    }
    else
    {
        r = sqrt(fma(a.x, a.x, a.y * a.y));
        a.x += r;
    //     absrz = 1 / sqrt((fma(r + a.x, r + a.x, a.y * a.y)) * r);
        absrz = 1 / sqrt((fma(a.x,a.x, a.y * a.y)) * r);
    //     res.x = absrz * (a.x + r);
        res.x = absrz * a.x;
        res.y = -absrz * a.y;
    }

    return res;
}


__inline__ __device__ double2 complexExponent(double2 a)
{
    double expr, sina, cosa;
    expr = exp(a.x);
    sincos(a.y, &sina, &cosa);
    a.x = expr * cosa;
    a.y = expr * sina;

    return a;
}

__inline__ __device__ double2 complexLog(double2 val)
{
    double2 logRes;
    logRes.x = log(sqrt(val.x * val.x + val.y * val.y));
    logRes.y = atan2(val.y, val.x);
    return logRes;
}

__inline__ __device__ double realComplexLog(double2 val)
{
    return log(sqrt(val.x * val.x + val.y * val.y));
}

__inline__ __device__ double2 cfma(double2 a, double2 b, double2 c)
{
    double2 res;

    res.x = fma(a.x,b.x,fma(-a.y,b.y,c.x));
    res.y = fma(a.x,b.y,fma(a.y,b.x,c.y));

    return res;
}

__inline__ __device__ double2 cfma(double a, double2 b, double2 c)
{
    double2 res;

    res.x = fma(a,b.x,c.x);
    res.y = fma(a,b.y,c.y);

    return res;
}

__inline__ __device__ double2 cfma(double a, double b, double2 c)
{
    double2 res;

    res.x = fma(a,b,c.x);
    res.y = c.y;

    return res;
}

__inline__ __device__ double2 cfmaconj(double2 a, double2 b, double2 c)
{
    double2 res;

    res.x = fma(a.x,b.x,fma(a.y,b.y,c.x));
    res.y = fma(a.x,-b.y,fma(a.y,b.x,c.y));

    return res;
}

__inline__ __device__ double logbesseli_05(double x)
{
    if(x == 0)
    {
        return 0.0;
    }
    if(x < 5)
    {
        return (log(SQRT_2_OVER_PI * sinh(x) / sqrt(x)));
    }
    return (x + log(SQRT_2_OVER_PI_OVER2 / sqrt(x)));
}

// FLOATS

__inline __device__ float2 buildComplex(float a, float b)
{
    float2 res;
    res.x = a; res.y = b;
    return res;
}

__inline __device__ float2 buildComplex(float a)
{
    return buildComplex(a,0.0f);
}

__inline __device__ float2 imagAdd(float2 a, float b)
{
    a.y += b;
    return a;
}

__inline __device__ float2 operator +(float2 a, float2 b)
{
    a.x += b.x;
    a.y += b.y;
    return a;
}

__inline __device__ float2 operator -(float2 a, float2 b)
{
    a.x -= b.x;
    a.y -= b.y;
    return a;
}

__inline __device__ float2 operator -(float a, float2 b)
{
    b.x = a - b.x;
    b.y = -b.y;
    return b;
}

__inline __device__ float2 operator +(float2 a, float b)
{
    a.x += b;
    return a;
}

__inline__ __device__ float2 operator /(float2 a, float2 b)
{
    float2 c;
    float denominator = 1 / (fmaf(b.x, b.x, b.y * b.y));
    c.x = (fmaf(a.x, b.x, a.y * b.y)) * denominator;
    c.y = (fmaf(a.y, b.x, -a.x * b.y)) * denominator;
        
//     c.x = (fma(a.x, b.x, a.y * b.y)) / (fma(b.x, b.x, b.y * b.y));
//     c.y = (fma(a.y, b.x, -a.x * b.y)) / (fma(b.x, b.x, b.y * b.y));

    return c;
}

__inline __device__ float2 operator /(float2 a, float b)
{
    a.x /= b;
    a.y /= b;
    return a;
}


__inline__ __device__ float2 operator /(float a, float2 b)
{
    a /= (fmaf(b.x, b.x, b.y * b.y));
    b.x *= a;
    b.y *= (-a);
    return b;
}


__inline__ __device__ float2 complexr(float2 b)
{
    float a = 1.0/(fmaf(b.x, b.x, b.y * b.y));
    b.x *= a;
    b.y *= (-a);
    return b;
}


__inline __device__ float2 operator *(float a, float2 b)
{
    b.x *= a;
    b.y *= a;
    return b;
}


__inline __device__ float2 operator *(float2 a, float2 b)
{
    float2 tmpVar = a;

    a.x = fmaf(a.x, b.x, -a.y * b.y);
    a.y = fmaf(tmpVar.x, b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ float realMult(float2 a, float2 b)
{
    return fmaf(a.x, b.x, -a.y * b.y);
}

__inline__  __device__ float2 conjMult(float2 a, float2 b)
{
    float2 tmpVar = a;

    a.x = fmaf(a.x, b.x, a.y * b.y);
    a.y = fmaf(tmpVar.x, -b.y, tmpVar.y * b.x);
    return a;
}

__inline __device__ float2 complexConj(float2 a)
{
    a.y *= -1.0f;
    return a;
}

__inline __device__ float2 complexSquare(float2 a)
{
    float2 res;
    res.x = (a.x - a.y) * (a.x + a.y);
//     res.x = fma(a.x,a.x,-a.y*a.y);
    res.y = 2.0f * a.x * a.y;
    return res;
}

__inline__ __device__ float2 complexSqrt(float2 a)
{
    float r, absrz;
    float2 res;

    if(a.y == 0.0f)
    {
        res.x = (a.x >= 0.0f ? sqrtf(a.x) : 0.0f);
        res.y = (a.x >= 0.0f ? 0.0f : sqrtf(-a.x));
    }
    else
    {
        r = sqrtf(fmaf(a.x, a.x, a.y * a.y));
        a.x = a.x + r;
    //     absrz = 1 / sqrt((fma(r + a.x, r + a.x, a.y * a.y)) / r);
        absrz = 1 / sqrtf((fmaf(a.x, a.x, a.y * a.y)) / r);
    //     res.x = absrz * (a.x + r);
        res.x = absrz * a.x;
        res.y = absrz * a.y;
    }
    return res;
}

__inline__ __device__ float2 rComplexSqrt(float2 a)
{
    float r, absrz;
    float2 res;
        
    if(a.y == 0.0f)
    {
        res.x = (a.x >= 0.0f ? 1/sqrtf(a.x) : 0.0f);
        res.y = (a.x >= 0.0f ? 0.0f : 1/sqrtf(-a.x));
    }
    else
    {
        r = sqrtf(fmaf(a.x, a.x, a.y * a.y));
        a.x += r;
    //     absrz = 1 / sqrt((fma(r + a.x, r + a.x, a.y * a.y)) * r);
        absrz = 1 / sqrtf((fmaf(a.x,a.x, a.y * a.y)) * r);
    //     res.x = absrz * (a.x + r);
        res.x = absrz * a.x;
        res.y = -absrz * a.y;
    }

    return res;
}


__inline__ __device__ float2 complexExponent(float2 a)
{
    float expr, sina, cosa;
    expr = expf(a.x);
    sincosf(a.y, &sina, &cosa);
    a.x = expr * cosa;
    a.y = expr * sina;

    return a;
}

__inline__ __device__ float2 complexLog(float2 val)
{
    float2 logRes;
    logRes.x = logf(sqrtf(val.x * val.x + val.y * val.y));
    logRes.y = atan2f(val.y, val.x);
    return logRes;
}

__inline__ __device__ float realComplexLog(float2 val)
{
    return log(sqrt(val.x * val.x + val.y * val.y));
}

__inline__ __device__ float2 cfma(float2 a, float2 b, float2 c)
{
    float2 res;

    res.x = fmaf(a.x,b.x,fmaf(-a.y,b.y,c.x));
    res.y = fmaf(a.x,b.y,fmaf(a.y,b.x,c.y));

    return res;
}

__inline__ __device__ float2 cfma(float a, float2 b, float2 c)
{
    float2 res;

    res.x = fmaf(a,b.x,c.x);
    res.y = fmaf(a,b.y,c.y);

    return res;
}

__inline__ __device__ float2 cfma(float a, float b, float2 c)
{
    float2 res;

    res.x = fmaf(a,b,c.x);
    res.y = c.y;

    return res;
}

__inline__ __device__ float2 cfmaconj(float2 a, float2 b, float2 c)
{
    float2 res;

    res.x = fmaf(a.x,b.x,fma(a.y,b.y,c.x));
    res.y = fmaf(a.x,-b.y,fma(a.y,b.x,c.y));

    return res;
}

__inline__ __device__ float logbesseli_05(float x)
{
    if(x == 0)
    {
        return 0.0f;
    }
    if(x < 5)
    {
        return (logf((float)(SQRT_2_OVER_PI) * sinhf(x) / sqrtf(x)));
    }
    return (x + logf((float)(SQRT_2_OVER_PI_OVER2) / sqrtf(x)));
}
