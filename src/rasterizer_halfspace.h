#pragma once

// Half-space / barycentric triangle rasterizer.
//
// Rasterizes in 2x2 pixel blocks, does perspective correct UV interpolation
// every 2nd column (e.g. at even X coordinates).
//
// References: Fabian Giesen's blog / pouet.net posts:
// - "The barycentric conspiracy" https://fgiesen.wordpress.com/2013/02/06/the-barycentric-conspirac/
// - "Triangle rasterization in practice" https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
// - "Optimizing the basic rasterizer" https://fgiesen.wordpress.com/2013/02/10/optimizing-the-basic-rasterizer/
// - "Fast software rasteriser in JavaScript?" https://www.pouet.net/topic.php?which=8760&page=1#c408170 and https://gist.github.com/rygorous/2486101
// - "Simple watertight triangle rasterizer" https://gist.github.com/rygorous/9b793cd21d876da928bf4c7f3e625908

#include <string.h>
#include "mathlib.h"
#include "platform.h"

#define SUBPIXEL_SHIFT (4)
#define SUBPIXEL_SCALE (1 << SUBPIXEL_SHIFT)

static inline int det2x2(int a, int b, int c, int d)
{
	int r = a * d - b * c;
	return r >> SUBPIXEL_SHIFT;
}

static inline int det2x2_fill_rule(int a, int b, int c, int d)
{
	int r = a * d + b * c; // Determinant
	if (c > 0 || (c == 0 && a > 0)) r++; // Top-left fill rule
	return (r - 1) >> SUBPIXEL_SHIFT;
}

static inline int to_fixed(float x)
{
	// -0.5f to place pixel centers at integer coords + 0.5
	// +0.5f afterwards is rounding factor
	return (int)((x - 0.5f) * SUBPIXEL_SCALE + 0.5f);
}

static int fixed_ceil(int x)
{
	return (x + SUBPIXEL_SCALE - 1) >> SUBPIXEL_SHIFT;
}

typedef struct raster_halfspace_t
{
	// bounding box in pixels
	int minx, miny, maxx, maxy;
	// edge vectors
	int dx12, dy12;
	int dx23, dy23;
	int dx31, dy31;
	// edge functions
	int c1, c2, c3;

	// 1.0/Z for each vertex
	float invz1, invz2, invz3;
	// UVs at each vertex, divided by Z
	float uv1x, uv1y;
	float uv2x, uv2y;
	float uv3x, uv3y;
} raster_halfspace_t;

typedef struct raster_halfspace_state_t
{
	int x, y;
	int c1, c2, c3;
	int cx1_00, cx2_00, cx3_00;

	bool first;
	float u_00, u_10, u_01, u_11, u_02, u_12;
	float v_00, v_10, v_01, v_11, v_02, v_12;

	float dudx, dvdx, dudy, dvdy;

	int cx1_10, cx2_10, cx3_10;
	int cx1_01, cx2_01, cx3_01;
	int cx1_11, cx2_11, cx3_11;
	bool inside_00, inside_01, inside_10, inside_11;
} raster_halfspace_state_t;


// Note: force inlining is needed on all these functions. Otherwise even with "static inline" the compiler
// decides to not inline some of them as soon as one than more call site is present, drastically reducing
// performance on Playdate.

FORCE_INLINE void raster_halfspace_begin(raster_halfspace_t* hs, const float3* p1, const float3* p2, const float3* p3, const float uvs[6])
{
	// convert coordinates to fixed point
	int x1 = to_fixed(p1->x), y1 = to_fixed(p1->y);
	int x2 = to_fixed(p2->x), y2 = to_fixed(p2->y);
	int x3 = to_fixed(p3->x), y3 = to_fixed(p3->y);

	// Note: do not check winding order or fully outside of screen; that is checked in calling code
	// check triangle winding order
	//int det = det2x2(x2 - x1, x3 - x1, y2 - y1, y3 - y1);
	//if (det == 0)
	//	return false; // zero area
	//if (det > 0)
	//	return false; // wrong winding

	// bounding box / clipping
	hs->minx = max2(fixed_ceil(min3(x1, x2, x3)), 0);
	hs->miny = max2(fixed_ceil(min3(y1, y2, y3)), 0);
	hs->maxx = min2(fixed_ceil(max3(x1, x2, x3)), SCREEN_X);
	hs->maxy = min2(fixed_ceil(max3(y1, y2, y3)), SCREEN_Y);
	//if (hs->minx >= hs->maxx || hs->miny >= hs->maxy)
	//	return false;

	// Edge vectors
	hs->dx12 = x1 - x2, hs->dy12 = y2 - y1;
	hs->dx23 = x2 - x3, hs->dy23 = y3 - y2;
	hs->dx31 = x3 - x1, hs->dy31 = y1 - y3;

	// We will rasterize 2x2 pixel blocks, so make sure starting coords are even
	hs->minx &= ~1;
	hs->miny &= ~1;

	// Edge functions
	int minx_fx = hs->minx << SUBPIXEL_SHIFT;
	int miny_fx = hs->miny << SUBPIXEL_SHIFT;
	hs->c1 = det2x2_fill_rule(hs->dx12, minx_fx - x1, hs->dy12, miny_fx - y1);
	hs->c2 = det2x2_fill_rule(hs->dx23, minx_fx - x2, hs->dy23, miny_fx - y2);
	hs->c3 = det2x2_fill_rule(hs->dx31, minx_fx - x3, hs->dy31, miny_fx - y3);

	hs->invz1 = 1.0f / p1->z, hs->invz2 = 1.0f / p2->z, hs->invz3 = 1.0f / p3->z;
	// Divide UVs by Z for perspective correction
	// https://www.scratchapixel.com/lessons/3d-basic-rendering/rasterization-practical-implementation/perspective-correct-interpolation-vertex-attributes.html
	hs->uv1x = uvs[0] * hs->invz1, hs->uv1y = uvs[1] * hs->invz1;
	hs->uv2x = uvs[2] * hs->invz2, hs->uv2y = uvs[3] * hs->invz2;
	hs->uv3x = uvs[4] * hs->invz3, hs->uv3y = uvs[5] * hs->invz3;
}

FORCE_INLINE void raster_halfspace_y_begin(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	memset(state, 0, sizeof(*state));

	state->y = hs->miny;
	state->c1 = hs->c1;
	state->c2 = hs->c2;
	state->c3 = hs->c3;
}

FORCE_INLINE bool raster_halfspace_y_continue(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	return state->y < hs->maxy;
}

FORCE_INLINE void raster_halfspace_y_step(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	state->c1 += hs->dx12 * 2;
	state->c2 += hs->dx23 * 2;
	state->c3 += hs->dx31 * 2;
	state->y += 2;
}

FORCE_INLINE void raster_halfspace_x_begin(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	state->cx1_00 = state->c1; state->cx2_00 = state->c2; state->cx3_00 = state->c3;
	state->first = true;
	state->x = hs->minx;
}

FORCE_INLINE bool raster_halfspace_x_continue(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	return state->x < hs->maxx;
}

FORCE_INLINE void raster_halfspace_x_step(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	state->cx1_00 += hs->dy12 * 2;
	state->cx2_00 += hs->dy23 * 2;
	state->cx3_00 += hs->dy31 * 2;

	state->x += 2;

	// Next loop iteration for x+0 UVs will use x+2 ones from current iteration
	state->u_00 = state->u_02;
	state->u_10 = state->u_12;
	state->v_00 = state->v_02;
	state->v_10 = state->v_12;
}

FORCE_INLINE int raster_halfspace_x_inner(const raster_halfspace_t* hs, raster_halfspace_state_t* state)
{
	state->cx1_10 = state->cx1_00 + hs->dx12, state->cx2_10 = state->cx2_00 + hs->dx23, state->cx3_10 = state->cx3_00 + hs->dx31;
	state->cx1_01 = state->cx1_00 + hs->dy12, state->cx2_01 = state->cx2_00 + hs->dy23, state->cx3_01 = state->cx3_00 + hs->dy31;
	state->cx1_11 = state->cx1_10 + hs->dy12, state->cx2_11 = state->cx2_10 + hs->dy23, state->cx3_11 = state->cx3_10 + hs->dy31;
	const int edges_00 = state->cx1_00 | state->cx2_00 | state->cx3_00;
	const int edges_01 = state->cx1_01 | state->cx2_01 | state->cx3_01;
	const int edges_10 = state->cx1_10 | state->cx2_10 | state->cx3_10;
	const int edges_11 = state->cx1_11 | state->cx2_11 | state->cx3_11;
	if ((edges_00 & edges_01 & edges_10 & edges_11) < 0) // all 4 pixels are outside
	{
		if (!state->first)
			return 1; // We already were inside triangle and got out, can skip the rest of the row
		else
			return 2; // Continue to next 2x2 block
	}

	state->inside_00 = edges_00 >= 0;
	state->inside_01 = edges_01 >= 0;
	state->inside_10 = edges_10 >= 0;
	state->inside_11 = edges_11 >= 0;

	// Barycentric coordinates for pixels.
	// Note: fast_rcp is same performance as division currently, with more artifacts. Keep division for now.
	// We only compute proper perspective correct barycentrics and UVs at each x2 horizontal step, i.e.
	// even pixels. The x+1 pixel interpolates between x+0 and x+2.
	if (state->first)
	{
		// First 2x2 block in this row, calculate UVs for x+0 pixels.
		const float bary_scale_00 = 1.0f / (state->cx1_00 * hs->invz3 + state->cx2_00 * hs->invz1 + state->cx3_00 * hs->invz2);
		const float bary_scale_10 = 1.0f / (state->cx1_10 * hs->invz3 + state->cx2_10 * hs->invz1 + state->cx3_10 * hs->invz2);
		const float bar_1_00 = state->cx1_00 * bary_scale_00;
		const float bar_1_10 = state->cx1_10 * bary_scale_10;
		const float bar_2_00 = state->cx2_00 * bary_scale_00;
		const float bar_2_10 = state->cx2_10 * bary_scale_10;
		const float bar_3_00 = state->cx3_00 * bary_scale_00;
		const float bar_3_10 = state->cx3_10 * bary_scale_10;
		state->u_00 = hs->uv3x * bar_1_00 + hs->uv1x * bar_2_00 + hs->uv2x * bar_3_00;
		state->u_10 = hs->uv3x * bar_1_10 + hs->uv1x * bar_2_10 + hs->uv2x * bar_3_10;
		state->v_00 = hs->uv3y * bar_1_00 + hs->uv1y * bar_2_00 + hs->uv2y * bar_3_00;
		state->v_10 = hs->uv3y * bar_1_10 + hs->uv1y * bar_2_10 + hs->uv2y * bar_3_10;
		state->first = false;
	}
	// Barycentrics and UVs for x+2 pixels.
	const int cx1_02 = state->cx1_01 + hs->dy12, cx2_02 = state->cx2_01 + hs->dy23, cx3_02 = state->cx3_01 + hs->dy31;
	const int cx1_12 = state->cx1_11 + hs->dy12, cx2_12 = state->cx2_11 + hs->dy23, cx3_12 = state->cx3_11 + hs->dy31;
	const float bary_scale_02 = 1.0f / (cx1_02 * hs->invz3 + cx2_02 * hs->invz1 + cx3_02 * hs->invz2);
	const float bary_scale_12 = 1.0f / (cx1_12 * hs->invz3 + cx2_12 * hs->invz1 + cx3_12 * hs->invz2);
	const float bar_1_02 = cx1_02 * bary_scale_02;
	const float bar_1_12 = cx1_12 * bary_scale_12;
	const float bar_2_02 = cx2_02 * bary_scale_02;
	const float bar_2_12 = cx2_12 * bary_scale_12;
	const float bar_3_02 = cx3_02 * bary_scale_02;
	const float bar_3_12 = cx3_12 * bary_scale_12;
	state->u_02 = hs->uv3x * bar_1_02 + hs->uv1x * bar_2_02 + hs->uv2x * bar_3_02;
	state->u_12 = hs->uv3x * bar_1_12 + hs->uv1x * bar_2_12 + hs->uv2x * bar_3_12;
	state->v_02 = hs->uv3y * bar_1_02 + hs->uv1y * bar_2_02 + hs->uv2y * bar_3_02;
	state->v_12 = hs->uv3y * bar_1_12 + hs->uv1y * bar_2_12 + hs->uv2y * bar_3_12;

	// Interpolate UVs for x+1 pixels.
	state->u_01 = (state->u_02 + state->u_00) * 0.5f;
	state->u_11 = (state->u_12 + state->u_10) * 0.5f;
	state->v_01 = (state->v_02 + state->v_00) * 0.5f;
	state->v_11 = (state->v_12 + state->v_10) * 0.5f;

	// UV derivatives
	state->dudx = state->u_01 - state->u_00;
	state->dvdx = state->v_01 - state->v_00;
	state->dudy = state->u_10 - state->u_00;
	state->dvdy = state->v_10 - state->v_00;
	return 0;
}
