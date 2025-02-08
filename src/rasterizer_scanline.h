#pragma once

// Scan-line triangle rasterizer.
//
// Based on Chris Hecker's "Perspective Texture Mapping" series (1995-1996)
// at https://chrishecker.com/Miscellaneous_Technical_Articles
//
// Implementation like in SUBAFXFL.CPP (subdividing affine mapper, fixed point
// edge step, float UV mapping).
//
// With additional tweaks:
// - Make it match raster rules of all the usual graphics APIs (Hecker's raster
//   rulesseem to be similar to D3D9).
// - Correctly handle spans that are clipped by edges of screen.

#include <assert.h>
#include <stdint.h>
#include "mathlib.h"
#include "platform.h"

typedef int32_t fixed28_4;
typedef int32_t fixed16_16;

FORCE_INLINE fixed28_4 FloatToFixed28_4(float v) { return (fixed28_4)(v * 16.0f); }
FORCE_INLINE float Fixed28_4ToFloat(fixed28_4 v) { return v * (1.0f / 16.0f); }
FORCE_INLINE fixed16_16 FloatToFixed16_16(float v) { return (fixed16_16)(v * 65536.0f); }
FORCE_INLINE float Fixed16_16ToFloat(fixed16_16 v) { return v * (1.0f / 65536.0f); }

static inline int Ceil28_4(fixed28_4 Value)
{
	int res;
	int numerator = Value - 1 + 16;
	if (numerator >= 0) {
		res = numerator / 16;
	}
	else {
		// deal with negative numerators correctly
		res = -((-numerator) / 16);
		res -= ((-numerator) % 16) ? 1 : 0;
	}
	return res;
}

static inline void FloorDivMod(int numerator, int denominator, int *r_floor, int *r_mod)
{
	assert(denominator > 0); // We assume it is positive
	if (numerator >= 0) {
		// positive case, C is okay
		*r_floor = numerator / denominator;
		*r_mod = numerator % denominator;
	}
	else {
		// numerator is negative, do the right thing
		*r_floor = -((-numerator) / denominator);
		*r_mod = (-numerator) % denominator;
		if (*r_mod) {
			// there is a remainder
			(*r_floor)--; *r_mod = denominator - *r_mod;
		}
	}
}

typedef struct tri_gradients {
	float invz[3];			// 1/z for each vertex
	float uz[3];			// u/z for each vertex
	float vz[3];			// v/z for each vertex
	float invz_dx, invz_dy;	// d(1/z)/dX, d(1/z)/dY
	float uz_dx, uz_dy;		// d(u/z)/dX, d(u/z)/dY
	float vz_dx, vz_dy;		// d(v/z)/dX, d(v/z)/dY
	float inv_det;
} tri_gradients;

static void tri_gradients_init(tri_gradients* t, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const float uv_scale)
{
	const float det = ((p1->x - p2->x) * (p0->y - p2->y)) - ((p0->x - p2->x) * (p1->y - p2->y));
	//if ((int)det == 0)
	//	return false; // zero area

	t->inv_det = 1.0f / det;
	float invdx = t->inv_det;
	float invdy = -t->inv_det;

	{
		float const invz = 1.0f / p0->z;
		t->invz[0] = invz;
		t->uz[0] = uvs[0] * (uv_scale * invz);
		t->vz[0] = uvs[1] * (uv_scale * invz);
	}
	{
		float const invz = 1.0f / p1->z;
		t->invz[1] = invz;
		t->uz[1] = uvs[2] * (uv_scale * invz);
		t->vz[1] = uvs[3] * (uv_scale * invz);
	}
	{
		float const invz = 1.0f / p2->z;
		t->invz[2] = invz;
		t->uz[2] = uvs[4] * (uv_scale * invz);
		t->vz[2] = uvs[5] * (uv_scale * invz);
	}

	float x02 = p0->x - p2->x;
	float y02 = p0->y - p2->y;
	float x12 = p1->x - p2->x;
	float y12 = p1->y - p2->y;

	t->invz_dx = invdx * (((t->invz[1] - t->invz[2]) * y02) - ((t->invz[0] - t->invz[2]) * y12));
	t->invz_dy = invdy * (((t->invz[1] - t->invz[2]) * x02) - ((t->invz[0] - t->invz[2]) * x12));

	t->uz_dx = invdx * (((t->uz[1] - t->uz[2]) * y02) - ((t->uz[0] - t->uz[2]) * y12));
	t->uz_dy = invdy * (((t->uz[1] - t->uz[2]) * x02) - ((t->uz[0] - t->uz[2]) * x12));

	t->vz_dx = invdx * (((t->vz[1] - t->vz[2]) * y02) - ((t->vz[0] - t->vz[2]) * y12));
	t->vz_dy = invdy * (((t->vz[1] - t->vz[2]) * x02) - ((t->vz[0] - t->vz[2]) * x12));
}

typedef struct tri_edge_uv {
	int X, x_step, numerator, denominator;	// DDA info for x
	int error_term;
	int Y, height;							// current y and vertical count
	float invz, invz_step, invz_step_extra;	// 1/z and step
	float uz, uz_step, uz_step_extra;		// u/z and step
	float vz, vz_step, vz_step_extra;		// v/z and step
	float x_prestep, y_prestep;
} tri_edge_uv;

typedef struct tri_edge {
	int X, x_step, numerator, denominator;	// DDA info for x
	int error_term;
	int Y, height;							// current y and vertical count
} tri_edge;

static void tri_edge_uv_init(tri_edge_uv* t, const tri_gradients* grad, const float3* p0, const float3* p1, const float3* p2, int top, int bottom)
{
	float3 vertTop = top == 0 ? *p0 : (top == 1 ? *p1 : *p2);
	float3 vertBottom = bottom == 0 ? *p0 : (bottom == 1 ? *p1 : *p2);
	fixed28_4 topXfx = FloatToFixed28_4(vertTop.x - 0.5f), botXfx = FloatToFixed28_4(vertBottom.x - 0.5f);
	fixed28_4 topYfx = FloatToFixed28_4(vertTop.y - 0.5f), botYfx = FloatToFixed28_4(vertBottom.y - 0.5f);
	t->Y = clamp_i(Ceil28_4(topYfx), 0, SCREEN_Y);
	const int YEnd = clamp_i(Ceil28_4(botYfx), 0, SCREEN_Y);
	t->height = YEnd - t->Y;
	assert(t->height >= 0);

	if (t->height)
	{
		fixed28_4 dN = botYfx - topYfx;
		fixed28_4 dM = botXfx - topXfx;

		fixed28_4 init_numerator = dM * 16 * t->Y - dM * topYfx + dN * topXfx - 1 + dN * 16;
		FloorDivMod(init_numerator, dN * 16, &t->X, &t->error_term);
		FloorDivMod(dM * 16, dN * 16, &t->x_step, &t->numerator);
		t->denominator = dN * 16;

		t->x_prestep = Fixed28_4ToFloat(t->X * 16 - topXfx);
		t->y_prestep = Fixed28_4ToFloat(t->Y * 16 - topYfx);

		t->invz = grad->invz[top] + t->y_prestep * grad->invz_dy + t->x_prestep * grad->invz_dx;
		t->invz_step = t->x_step * grad->invz_dx + grad->invz_dy;
		t->invz_step_extra = grad->invz_dx;

		t->uz = grad->uz[top] + t->y_prestep * grad->uz_dy + t->x_prestep * grad->uz_dx;
		t->uz_step = t->x_step * grad->uz_dx + grad->uz_dy;
		t->uz_step_extra = grad->uz_dx;

		t->vz = grad->vz[top] + t->y_prestep * grad->vz_dy + t->x_prestep * grad->vz_dx;
		t->vz_step = t->x_step * grad->vz_dx + grad->vz_dy;
		t->vz_step_extra = grad->vz_dx;
	}
}

static void tri_edge_init(tri_edge* t, const float3* p0, const float3* p1, const float3* p2, int top, int bottom)
{
	float3 vertTop = top == 0 ? *p0 : (top == 1 ? *p1 : *p2);
	float3 vertBottom = bottom == 0 ? *p0 : (bottom == 1 ? *p1 : *p2);
	fixed28_4 topXfx = FloatToFixed28_4(vertTop.x - 0.5f), botXfx = FloatToFixed28_4(vertBottom.x - 0.5f);
	fixed28_4 topYfx = FloatToFixed28_4(vertTop.y - 0.5f), botYfx = FloatToFixed28_4(vertBottom.y - 0.5f);
	t->Y = clamp_i(Ceil28_4(topYfx), 0, SCREEN_Y);
	const int YEnd = clamp_i(Ceil28_4(botYfx), 0, SCREEN_Y);
	t->height = YEnd - t->Y;
	assert(t->height >= 0);

	if (t->height)
	{
		fixed28_4 dN = botYfx - topYfx;
		fixed28_4 dM = botXfx - topXfx;

		fixed28_4 init_numerator = dM * 16 * t->Y - dM * topYfx + dN * topXfx - 1 + dN * 16;
		FloorDivMod(init_numerator, dN * 16, &t->X, &t->error_term);
		FloorDivMod(dM * 16, dN * 16, &t->x_step, &t->numerator);
		t->denominator = dN * 16;
	}
}

static inline void tri_edge_uv_step(tri_edge_uv* t)
{
	t->X += t->x_step; t->Y++; t->height--;
	t->uz += t->uz_step; t->vz += t->vz_step; t->invz += t->invz_step;

	t->error_term += t->numerator;
	if (t->error_term >= t->denominator) {
		t->X++;
		t->error_term -= t->denominator;
		t->invz += t->invz_step_extra;
		t->uz += t->uz_step_extra; t->vz += t->vz_step_extra;
	}
}

static inline void tri_edge_step(tri_edge* t)
{
	t->X += t->x_step; t->Y++; t->height--;

	t->error_term += t->numerator;
	if (t->error_term >= t->denominator) {
		t->X++;
		t->error_term -= t->denominator;
	}
}

#define kAffineLength (8)

typedef struct raster_scanline_line_t {
	float invz_left, uz_left, vz_left;
	float invz_aff_step, uz_aff_step ,vz_aff_step;
	float invz_right, uz_right, vz_right;

	float z_left, u_left, v_left;
	float z_right, u_right, v_right;

	fixed16_16 U, V, dU, dV;

	int aff_spans, aff_remainder;

	float clip_adj_l, clip_adj_r;
	int x_start, x_end;
} raster_scanline_line_t;

FORCE_INLINE bool raster_scanline_line_init(raster_scanline_line_t* line, const tri_gradients* grad, const tri_edge_uv* edge_l, const tri_edge_uv* edge_r)
{
	assert(edge_l->Y >= 0);
	assert(edge_l->Y < SCREEN_Y);
	assert(edge_l->Y == edge_r->Y);

	line->x_start = max2(0, edge_l->X);
	line->x_end = min2(SCREEN_X, edge_r->X);
	if (line->x_end <= line->x_start)
		return false;

	line->clip_adj_l = (float)(line->x_start - edge_l->X);
	line->clip_adj_r = (float)(edge_r->X - line->x_end + 1);

	line->invz_left = edge_l->invz + line->clip_adj_l * grad->invz_dx;
	line->uz_left = edge_l->uz + line->clip_adj_l * grad->uz_dx;
	line->vz_left = edge_l->vz + line->clip_adj_l * grad->vz_dx;

	line->invz_aff_step = grad->invz_dx * kAffineLength;
	line->uz_aff_step = grad->uz_dx * kAffineLength;
	line->vz_aff_step = grad->vz_dx * kAffineLength;

	line->invz_right = line->invz_left + line->invz_aff_step;
	line->uz_right = line->uz_left + line->uz_aff_step;
	line->vz_right = line->vz_left + line->vz_aff_step;

	line->z_left = 1.0f / line->invz_left;
	line->u_left = line->z_left * line->uz_left;
	line->v_left = line->z_left * line->vz_left;

	const int width = line->x_end - line->x_start;
	line->aff_spans = width / kAffineLength;
	line->aff_remainder = width % kAffineLength;
	if (!line->aff_remainder)
	{
		line->aff_spans--;
		line->aff_remainder = kAffineLength;
	}

	line->U = line->V = line->dU = line->dV = 0;
	return true;
}

FORCE_INLINE bool raster_scanline_spans_continue(raster_scanline_line_t* line)
{
	return line->aff_spans-- > 0;
}

FORCE_INLINE void raster_scanline_spans_step(raster_scanline_line_t* line)
{
	line->z_left = line->z_right;
	line->u_left = line->u_right;
	line->v_left = line->v_right;

	line->invz_right += line->invz_aff_step;
	line->uz_right += line->uz_aff_step;
	line->vz_right += line->vz_aff_step;
}

FORCE_INLINE void raster_scanline_inner_begin(raster_scanline_line_t* line)
{
	line->z_right = 1.0f / line->invz_right;
	line->u_right = line->z_right * line->uz_right;
	line->v_right = line->z_right * line->vz_right;

	line->U = FloatToFixed16_16(line->u_left);
	line->V = FloatToFixed16_16(line->v_left);
	line->dU = FloatToFixed16_16(line->u_right - line->u_left) / kAffineLength;
	line->dV = FloatToFixed16_16(line->v_right - line->v_left) / kAffineLength;
}

FORCE_INLINE void raster_scanline_inner_step(raster_scanline_line_t* line)
{
	line->U += line->dU;
	line->V += line->dV;
}

FORCE_INLINE bool raster_scanline_has_rem(raster_scanline_line_t* line)
{
	return line->aff_remainder;
}

FORCE_INLINE void raster_scanline_rem_begin(raster_scanline_line_t* line, const tri_gradients* grad, const tri_edge_uv* edge_r)
{
	line->z_right = 1.0f / (edge_r->invz - grad->invz_dx * line->clip_adj_r);
	line->u_right = line->z_right * (edge_r->uz - grad->uz_dx * line->clip_adj_r);
	line->v_right = line->z_right * (edge_r->vz - grad->vz_dx * line->clip_adj_r);

	line->U = FloatToFixed16_16(line->u_left);
	line->V = FloatToFixed16_16(line->v_left);
	if (--line->aff_remainder)
	{
		// guard against div-by-0 for 1 pixel lines
		line->dU = FloatToFixed16_16(line->u_right - line->u_left) / line->aff_remainder;
		line->dV = FloatToFixed16_16(line->v_right - line->v_left) / line->aff_remainder;
	}
}

static bool tri_sort_vertices(const float3* p0, const float3* p1, const float3* p2, int* r_top, int* r_mid, int* r_bottom)
{
	int v_mid_cmp, v_bot_cmp;
	float Y0 = p0->y, Y1 = p1->y, Y2 = p2->y;

	// sort vertices in y
	if (Y0 < Y1) {
		if (Y2 < Y0) {
			*r_top = 2; *r_mid = 0; *r_bottom = 1;
			v_mid_cmp = 0; v_bot_cmp = 1;
		}
		else {
			*r_top = 0;
			if (Y1 < Y2) {
				*r_mid = 1; *r_bottom = 2;
				v_mid_cmp = 1; v_bot_cmp = 2;
			}
			else {
				*r_mid = 2; *r_bottom = 1;
				v_mid_cmp = 2; v_bot_cmp = 1;
			}
		}
	}
	else {
		if (Y2 < Y1) {
			*r_top = 2; *r_mid = 1; *r_bottom = 0;
			v_mid_cmp = 1; v_bot_cmp = 0;
		}
		else {
			*r_top = 1;
			if (Y0 < Y2) {
				*r_mid = 0; *r_bottom = 2;
				v_mid_cmp = 3; v_bot_cmp = 2;
			}
			else {
				*r_mid = 2; *r_bottom = 0;
				v_mid_cmp = 2; v_bot_cmp = 3;
			}
		}
	}
	// the triangle is anti-clockwise, so if bottom > middle then middle is right
	return v_bot_cmp <= v_mid_cmp;
}

typedef struct raster_scanline_simple_t
{
	tri_edge e_top_bottom, e_top_mid, e_mid_bottom;
	tri_edge* e_left;
	tri_edge* e_right;
	bool mid_is_left;
} raster_scanline_simple_t;


typedef struct raster_scanline_uv_t
{
	tri_gradients grad;
	tri_edge_uv e_top_bottom, e_top_mid, e_mid_bottom;
	tri_edge_uv* e_left;
	tri_edge_uv* e_right;
	bool mid_is_left;
	int v_top, v_mid, v_bottom;
} raster_scanline_uv_t;

FORCE_INLINE void raster_scanline_simple_begin(raster_scanline_simple_t* st, const float3* p0, const float3* p1, const float3* p2)
{
	int v_top, v_mid, v_bottom;
	st->mid_is_left = tri_sort_vertices(p0, p1, p2, &v_top, &v_mid, &v_bottom);

	// triangle outside of Y range?
	//const float y_top = v_top == 0 ? p0->y : (v_top == 1 ? p1->y : p2->y);
	//if (y_top >= SCREEN_Y) return false;
	//const float y_bottom = v_bottom == 0 ? p0->y : (v_bottom == 1 ? p1->y : p2->y);
	//if (y_bottom < 0.0f) return false;

	// zero area?
	//const float det = ((p1->x - p2->x) * (p0->y - p2->y)) - ((p0->x - p2->x) * (p1->y - p2->y));
	//if ((int)det == 0)
	//	return false; 

	tri_edge_init(&st->e_top_bottom, p0, p1, p2, v_top, v_bottom);
	tri_edge_init(&st->e_top_mid, p0, p1, p2, v_top, v_mid);
	tri_edge_init(&st->e_mid_bottom, p0, p1, p2, v_mid, v_bottom);

	st->e_left = st->e_right = NULL;
}

FORCE_INLINE void raster_scanline_uv_begin(raster_scanline_uv_t* st, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const float uv_scale)
{
	st->mid_is_left = tri_sort_vertices(p0, p1, p2, &st->v_top, &st->v_mid, &st->v_bottom);

	// triangle outside of Y range?
	//const float y_top = st->v_top == 0 ? p0->y : (st->v_top == 1 ? p1->y : p2->y);
	//if (y_top >= SCREEN_Y) return false;
	//const float y_bottom = st->v_bottom == 0 ? p0->y : (st->v_bottom == 1 ? p1->y : p2->y);
	//if (y_bottom < 0.0f) return false;

	tri_gradients_init(&st->grad, p0, p1, p2, uvs, uv_scale);
	tri_edge_uv_init(&st->e_top_bottom, &st->grad, p0, p1, p2, st->v_top, st->v_bottom);
	tri_edge_uv_init(&st->e_top_mid, &st->grad, p0, p1, p2, st->v_top, st->v_mid);
	tri_edge_uv_init(&st->e_mid_bottom, &st->grad, p0, p1, p2, st->v_mid, st->v_bottom);

	st->e_left = st->e_right = NULL;
}

FORCE_INLINE int raster_scanline_simple_begin1(raster_scanline_simple_t* st)
{
	st->e_left = st->mid_is_left ? &st->e_top_bottom : &st->e_top_mid;
	st->e_right = st->mid_is_left ? &st->e_top_mid : &st->e_top_bottom;
	return st->e_top_mid.height;
}

FORCE_INLINE int raster_scanline_uv_begin1(raster_scanline_uv_t* st)
{
	st->e_left = st->mid_is_left ? &st->e_top_bottom : &st->e_top_mid;
	st->e_right = st->mid_is_left ? &st->e_top_mid : &st->e_top_bottom;
	return st->e_top_mid.height;
}

FORCE_INLINE int raster_scanline_simple_begin2(raster_scanline_simple_t* st)
{
	st->e_left = st->mid_is_left ? &st->e_top_bottom : &st->e_mid_bottom;
	st->e_right = st->mid_is_left ? &st->e_mid_bottom : &st->e_top_bottom;
	return st->e_mid_bottom.height;
}

FORCE_INLINE int raster_scanline_uv_begin2(raster_scanline_uv_t* st)
{
	st->e_left = st->mid_is_left ? &st->e_top_bottom : &st->e_mid_bottom;
	st->e_right = st->mid_is_left ? &st->e_mid_bottom : &st->e_top_bottom;
	return st->e_mid_bottom.height;
}

FORCE_INLINE void raster_scanline_simple_y_step(raster_scanline_simple_t* st)
{
	tri_edge_step(st->e_left); tri_edge_step(st->e_right);
}

FORCE_INLINE void raster_scanline_uv_y_step(raster_scanline_uv_t* st)
{
	tri_edge_uv_step(st->e_left); tri_edge_uv_step(st->e_right);
}
