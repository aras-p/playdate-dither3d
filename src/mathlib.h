// SPDX-License-Identifier: Unlicense

#pragma once

#include <stdint.h>
#include <stdbool.h>
#include <math.h>


#if defined(_MSC_VER)
#  define FORCE_INLINE static __forceinline
#else
#  define FORCE_INLINE static inline __attribute__((always_inline)) __attribute__((__unused__))
#endif


#if TARGET_PLAYDATE
static inline float _lib3d_sqrtf(float x)
{
	float result;
	asm ("vsqrt.f32 %0, %1" : "=t"(result) : "t"(x));
	return result;
}
#define sqrtf(f) _lib3d_sqrtf(f)
#endif

static inline uint32_t swap(uint32_t n)
{
#if TARGET_PLAYDATE
	//return __REV(n);
	uint32_t result;
	__asm volatile ("rev %0, %1" : "=l" (result) : "l" (n));
	return(result);
#else
	return ((n & 0xff000000) >> 24) | ((n & 0xff0000) >> 8) | ((n & 0xff00) << 8) | (n << 24);
#endif
}

static inline int highest_bit(uint32_t x)
{
	if (x == 0)
		return -1;
#if _MSC_VER
	unsigned long res = 0;
	_BitScanReverse(&res, x);
	return (int)res;
#elif __clang__ || __GNUC__ || __GCC__
	return 31 - __builtin_clz(x);
#else
#error Unknown platform
#endif
}

// https://en.wikipedia.org/wiki/Fast_inverse_square_root
inline float fast_rsqrt(float x)
{
	union {
		float f;
		uint32_t u;
	} fu;
	fu.f = x;
	fu.u = 0x5f3759df - (fu.u >> 1);
	fu.f *= 1.5f - (x * 0.5f * fu.f * fu.f);
	return fu.f;
}
inline float fast_rcp(float x) { x = fast_rsqrt(x); return x * x; }
inline float fast_div(float x, float y) { return x * fast_rsqrt(y); }
inline float fast_sqrt(float x) { return x * fast_rsqrt(x); }

static inline int min2(int a, int b)
{
	return a < b ? a : b;
}

static inline int max2(int a, int b)
{
	return a > b ? a : b;
}

static inline float min2f(float a, float b)
{
	return a < b ? a : b;
}

static inline float max2f(float a, float b)
{
	return a > b ? a : b;
}

static inline int min3(int a, int b, int c)
{
	return (a < b) ? ((a < c) ? a : c) : ((b < c) ? b : c);
}

static inline int max3(int a, int b, int c)
{
	return (a > b) ? ((a > c) ? a : c) : ((b > c) ? b : c);
}

static inline int clamp_i(int v, int a, int b)
{
	if (v < a) v = a;
	if (v > b) v = b;
	return v;
}


#define M_PIf (3.14159265f)

static inline float saturate(float v)
{
	if (v < 0.0f)
		v = 0.0f;
	if (v > 1.0f)
		v = 1.0f;
	return v;
}

static inline float fract(float v)
{
	return v - floorf(v);
}

// https://gist.github.com/volkansalma/2972237?permalink_comment_id=3872525#gistcomment-3872525
static inline float atan2f_approx(float y, float x)
{
	float absy = fabsf(y) + 1.0e-8f; // avoid 0/0
	float r = (x - copysignf(absy, x)) / (absy + fabsf(x));
	float a = M_PIf * 0.5f - copysignf(M_PIf * 0.25f, x);
	a += (0.1963f * r * r - 0.9817f) * r;
	return copysignf(a, y);
}

static inline float lerp(float a, float b, float t)
{
	return a * (1.0f - t) + b * t;
}

static inline float smoothstep(float a, float b, float v)
{
	float t = (v - a) / (b - a);
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;
	return t * t * (3.0f - 2.0f * t);
}

static inline float invlerp(float a, float b, float v)
{
	float t = (v - a) / (b - a);
	if (t < 0.0f) t = 0.0f;
	if (t > 1.0f) t = 1.0f;
	return t;
}

static inline uint32_t XorShift32(uint32_t* state)
{
	uint32_t x = *state;
	x ^= x << 13;
	x ^= x >> 17;
	x ^= x << 15;
	*state = x;
	return x;
}

static inline float RandomFloat01(uint32_t* state)
{
	return (XorShift32(state) & 0xFFFFFF) / 16777216.0f;
}

typedef struct float2
{
	float x;
	float y;
} float2;

typedef struct float3
{
	float x;
	float y;
	float z;
} float3;

static inline float3 f3(float x, float y, float z)
{
	return (float3) {.x=x, .y=y, .z=z };
}

static inline float3 v3_cross(float3 a, float3 b)
{
	return (float3){ .x = a.y * b.z - a.z * b.y, .y = a.z * b.x - a.x * b.z, .z = a.x * b.y - a.y * b.x };
}

static inline float v3_dot(float3 a, float3 b)
{
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

static inline float3 v3_normalize(float3 v)
{
	float d = sqrtf(v.x * v.x + v.y * v.y + v.z * v.z);
	float id = 1.0f / d;
	float3 r = { v.x * id, v.y * id, v.z * id };
	return r;
}


static inline float v3_lensq(float3* v)
{
	return v->x * v->x + v->y * v->y + v->z * v->z;
}

static inline float v3_len(float3* v)
{
	return sqrtf(v3_lensq(v));
}

static inline float3 v3_add(float3 a, float3 v)
{
	return (float3) { a.x + v.x, a.y + v.y, a.z + v.z };
}
static inline float3 v3_sub(float3 a, float3 v)
{
	return (float3) { a.x - v.x, a.y - v.y, a.z - v.z };
}
static inline float3 v3_subfl(float3 a, float v)
{
	return (float3) { a.x - v, a.y - v, a.z - v };
}
static inline float3 v3_mul(float3 a, float3 b)
{
	return (float3) { a.x * b.x, a.y * b.y, a.z * b.z };
}
static inline float3 v3_mulfl(float3 v, float s)
{
	return (float3) { v.x*s, v.y*s, v.z*s };
}
static inline float3 v3_lerp(float3 a, float3 b, float f)
{
	float ff = 1.0f - f;
	float3 r = { a.x*ff + b.x*f, a.y*ff + b.y*f, a.z*ff + b.z*f };
	return r;
}

static inline float3 v3_fract(float3 v)
{
	float3 r = { fract(v.x), fract(v.y), fract(v.z) };
	return r;
}

static inline float3 v3_abs(float3 v)
{
	float3 r = { fabsf(v.x), fabsf(v.y), fabsf(v.z) };
	return r;
}

static inline float3 v3_reflect(float3 v, float3 n)
{
	return v3_sub(v, v3_mulfl(n, 2 * v3_dot(v, n)));
}

float3 v3_tri_normal(const float3* p1, const float3* p2, const float3* p3);

typedef struct xform
{
	float x;
	float y;
	float z;
	float m[3][3];
} xform;
extern xform xform_identity;
xform xform_make(float m11, float m12, float m13, float m21, float m22, float m23, float m31, float m32, float m33);
xform xform_axes(float3 ax, float3 ay, float3 az, float3 pos);
xform xform_make_translate(float dx, float dy, float dz);
xform xform_make_axis_angle(float angle, float3 axis);
xform xform_multiply(const xform* l, const xform* r);
float3 xform_transform_pt(const xform* m, float3 p);
float xform_get_determinant(xform* m);
