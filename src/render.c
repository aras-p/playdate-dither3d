#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "render.h"
#include "mathlib.h"
#include "platform.h"
#include "util/pixel_ops.h"

// Parts of the scanline rasterizer are originally based on:
// - mini3d from Playdate SDK examples. Created by Dave Hayden on 10/20/15. Copyright © 2015 Panic, Inc. All rights reserved.


FORCE_INLINE void _drawMaskPattern(uint32_t* p, uint32_t mask, uint32_t color)
{
	if ( mask == 0xffffffff )
		*p = color;
	else
		*p = (*p & ~mask) | (color & mask);
}

static void
drawFragment(uint32_t* row, int x1, int x2, uint32_t color)
{
	if ( x2 < 0 || x1 >= SCREEN_X )
		return;
	
	if ( x1 < 0 )
		x1 = 0;
	
	if ( x2 > SCREEN_X )
		x2 = SCREEN_X;
	
	if ( x1 > x2 )
		return;
	
	// Operate on 32 bits at a time
	
	int startbit = x1 % 32;
	uint32_t startmask = swap((1 << (32 - startbit)) - 1);
	int endbit = x2 % 32;
	uint32_t endmask = swap(((1 << endbit) - 1) << (32 - endbit));
	
	int col = x1 / 32;
	uint32_t* p = row + col;

	if ( col == x2 / 32 )
	{
		uint32_t mask = 0;
		
		if ( startbit > 0 && endbit > 0 )
			mask = startmask & endmask;
		else if ( startbit > 0 )
			mask = startmask;
		else if ( endbit > 0 )
			mask = endmask;
		
		_drawMaskPattern(p, mask, color);
	}
	else
	{
		int x = x1;
		
		if ( startbit > 0 )
		{
			_drawMaskPattern(p++, startmask, color);
			x += (32 - startbit);
		}
		
		while ( x + 32 <= x2 )
		{
			_drawMaskPattern(p++, 0xffffffff, color);
			x += 32;
		}
		
		if ( endbit > 0 )
			_drawMaskPattern(p, endmask, color);
	}
}

static inline int32_t slope(float x1, float y1, float x2, float y2, int shift)
{
	float dx = x2-x1;
	float dy = y2-y1;
	
	if ( dy < 1 )
		return (int32_t)(dx * (1<< shift));
	else
		return (int32_t)(dx / dy * (1<<shift));
}

void drawLine(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, int thick, const uint8_t pattern[8])
{
	if ( p1->y > p2->y )
	{
		const float3* tmp = p1;
		p1 = p2;
		p2 = tmp;
	}

	int y = (int)p1->y;
	int endy = (int)p2->y;
	
	if ( y >= SCREEN_Y || endy < 0 || min2f(p1->x, p2->x) >= SCREEN_X || max2f(p1->x, p2->x) < 0 )
		return;
	
	int32_t x = (int32_t)(p1->x * (1<<16));
	int32_t dx = slope(p1->x, p1->y, p2->x, p2->y, 16);
	float py = p1->y;
	
	if ( y < 0 )
	{
		x += (int32_t)(-p1->y * dx);
		y = 0;
		py = 0;
	}

	int32_t x1 = (int32_t)(x + dx * (y+1-py));

	while ( y <= endy )
	{
		uint8_t p = pattern[y%8];
		uint32_t color = (p<<24) | (p<<16) | (p<<8) | p;
		
		if ( y == endy )
			x1 = (int32_t)(p2->x * (1<<16));
		
		if ( dx < 0 )
			drawFragment((uint32_t*)&bitmap[y*rowstride], x1>>16, (x>>16) + thick, color);
		else
			drawFragment((uint32_t*)&bitmap[y*rowstride], x>>16, (x1>>16) + thick, color);
		
		if ( ++y == SCREEN_Y )
			break;

		x = x1;
		x1 += dx;
	}
}

// --------------------------------------------------------------------------
// Based on Chris Hecker's "Perspective Texture Mapping" series (1995-1996)
// at https://chrishecker.com/Miscellaneous_Technical_Articles
// Implementation like in SUBAFXFL.CPP (subdividing affine mapper, fixed point
// edge step, float UV mapping). With additional tweaks:
// - Make it match raster rules of all the usual graphics APIs (Hecker's raster
//   rulesseem to be similar to D3D9).
// - Correctly handle spans that are clipped by edges of screen.

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

static bool tri_gradients_init(tri_gradients* t, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const float uv_scale)
{
	const float det = ((p1->x - p2->x) * (p0->y - p2->y)) - ((p0->x - p2->x) * (p1->y - p2->y));
	if ((int)det == 0)
		return false; // zero area

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
	return true;
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

FORCE_INLINE bool raster_scanline_simple_begin(raster_scanline_simple_t* st, const float3* p0, const float3* p1, const float3* p2)
{
	int v_top, v_mid, v_bottom;
	st->mid_is_left = tri_sort_vertices(p0, p1, p2, &v_top, &v_mid, &v_bottom);

	// triangle outside of Y range?
	const float y_top = v_top == 0 ? p0->y : (v_top == 1 ? p1->y : p2->y);
	if (y_top >= SCREEN_Y) return false;
	const float y_bottom = v_bottom == 0 ? p0->y : (v_bottom == 1 ? p1->y : p2->y);
	if (y_bottom < 0.0f) return false;

	// zero area?
	const float det = ((p1->x - p2->x) * (p0->y - p2->y)) - ((p0->x - p2->x) * (p1->y - p2->y));
	if ((int)det == 0)
		return false; 

	tri_edge_init(&st->e_top_bottom, p0, p1, p2, v_top, v_bottom);
	tri_edge_init(&st->e_top_mid, p0, p1, p2, v_top, v_mid);
	tri_edge_init(&st->e_mid_bottom, p0, p1, p2, v_mid, v_bottom);

	st->e_left = st->e_right = NULL;
	return true;
}

FORCE_INLINE bool raster_scanline_uv_begin(raster_scanline_uv_t* st, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const float uv_scale)
{
	st->mid_is_left = tri_sort_vertices(p0, p1, p2, &st->v_top, &st->v_mid, &st->v_bottom);

	// triangle outside of Y range?
	const float y_top = st->v_top == 0 ? p0->y : (st->v_top == 1 ? p1->y : p2->y);
	if (y_top >= SCREEN_Y) return false;
	const float y_bottom = st->v_bottom == 0 ? p0->y : (st->v_bottom == 1 ? p1->y : p2->y);
	if (y_bottom < 0.0f) return false;

	if (!tri_gradients_init(&st->grad, p0, p1, p2, uvs, uv_scale))
		return false; // zero area
	tri_edge_uv_init(&st->e_top_bottom, &st->grad, p0, p1, p2, st->v_top, st->v_bottom);
	tri_edge_uv_init(&st->e_top_mid, &st->grad, p0, p1, p2, st->v_top, st->v_mid);
	tri_edge_uv_init(&st->e_mid_bottom, &st->grad, p0, p1, p2, st->v_mid, st->v_bottom);

	st->e_left = st->e_right = NULL;
	return true;
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

// --------------------------------------------------------------------------

FORCE_INLINE void draw_scanline_pattern(uint8_t* bitmap, int rowstride, const tri_edge* edge_l, const tri_edge* edge_r, const uint8_t* pattern)
{
	int x = max2(0, edge_l->X);
	int endx = min2(SCREEN_X, edge_r->X);
	if (endx <= x) return;

	uint8_t* row = bitmap + edge_l->Y * rowstride;

	const uint8_t pat = pattern[edge_l->Y % 8];
	uint32_t pattern4 = (pat << 24) | (pat << 16) | (pat << 8) | pat;

	// write out pixels 32 at a time
	int startbit = x % 32;
	uint32_t startmask = swap((1 << (32 - startbit)) - 1);
	int endbit = endx % 32;
	uint32_t endmask = swap(((1 << endbit) - 1) << (32 - endbit));

	int col = x / 32;
	uint32_t* p = (uint32_t*)row + col;

	if (col == endx / 32)
	{
		uint32_t mask = 0;
		if (startbit > 0 && endbit > 0)
			mask = startmask & endmask;
		else if (startbit > 0)
			mask = startmask;
		else if (endbit > 0)
			mask = endmask;
		_drawMaskPattern(p, mask, pattern4);
	}
	else
	{
		if (startbit > 0)
		{
			_drawMaskPattern(p++, startmask, pattern4);
			x += (32 - startbit);
		}
		while (x + 32 <= endx)
		{
			_drawMaskPattern(p++, 0xffffffff, pattern4);
			x += 32;
		}
		if (endbit > 0)
		{
			_drawMaskPattern(p, endmask, pattern4);
		}
	}
}

static void draw_triangle_pattern_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t* pattern)
{
	raster_scanline_simple_t tri;
	if (!raster_scanline_simple_begin(&tri, p0, p1, p2))
		return;

	int height = raster_scanline_simple_begin1(&tri);
	while (height--) {
		draw_scanline_pattern(bitmap, rowstride, tri.e_left, tri.e_right, pattern);
		raster_scanline_simple_y_step(&tri);
	}
	height = raster_scanline_simple_begin2(&tri);
	while (height--) {
		draw_scanline_pattern(bitmap, rowstride, tri.e_left, tri.e_right, pattern);
		raster_scanline_simple_y_step(&tri);
	}
}

FORCE_INLINE void draw_scanline_bluenoise(uint8_t* bitmap, int rowstride, const tri_edge* edge_l, const tri_edge* edge_r, const uint8_t tri_color)
{
	int x = max2(0, edge_l->X);
	int endx = min2(SCREEN_X, edge_r->X);
	if (endx <= x) return;

	const uint8_t* noise_row = get_blue_noise_buffer() + edge_l->Y * SCREEN_X;
	uint8_t* row = bitmap + edge_l->Y * rowstride;

	// write out pixels 32 at a time
	uint32_t mask = 0;
	uint32_t* p = (uint32_t*)row + x / 32;
	uint32_t color = 0;

	while (x < endx)
	{
		mask |= 0x80000000u >> (x & 31);
		uint32_t pix = tri_color > noise_row[x] ? 0x80000000 : 0;
		color |= pix >> (x & 31);
		x++;
		if (x % 32 == 0)
		{
			_drawMaskPattern(p++, swap(mask), swap(color));
			mask = 0;
			color = 0;
		}
	}
	_drawMaskPattern(p, swap(mask), swap(color));
}

static void draw_triangle_bluenoise_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t tri_color)
{
	raster_scanline_simple_t tri;
	if (!raster_scanline_simple_begin(&tri, p0, p1, p2))
		return;

	int height = raster_scanline_simple_begin1(&tri);
	while (height--) {
		draw_scanline_bluenoise(bitmap, rowstride, tri.e_left, tri.e_right, tri_color);
		raster_scanline_simple_y_step(&tri);
	}
	height = raster_scanline_simple_begin2(&tri);
	while (height--) {
		draw_scanline_bluenoise(bitmap, rowstride, tri.e_left, tri.e_right, tri_color);
		raster_scanline_simple_y_step(&tri);
	}
}


// --------------------------------------------------------------------------
// Half-space / barycentric triangle rasterizer

// References: Fabian Giesen's blog / pouet.net posts:
// - "The barycentric conspiracy" https://fgiesen.wordpress.com/2013/02/06/the-barycentric-conspirac/
// - "Triangle rasterization in practice" https://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
// - "Optimizing the basic rasterizer" https://fgiesen.wordpress.com/2013/02/10/optimizing-the-basic-rasterizer/
// - "Fast software rasteriser in JavaScript?" https://www.pouet.net/topic.php?which=8760&page=1#c408170 and https://gist.github.com/rygorous/2486101
// - "Simple watertight triangle rasterizer" https://gist.github.com/rygorous/9b793cd21d876da928bf4c7f3e625908

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

FORCE_INLINE bool raster_halfspace_begin(raster_halfspace_t* hs, const float3* p1, const float3* p2, const float3* p3, const float uvs[6])
{
	// convert coordinates to fixed point
	int x1 = to_fixed(p1->x), y1 = to_fixed(p1->y);
	int x2 = to_fixed(p2->x), y2 = to_fixed(p2->y);
	int x3 = to_fixed(p3->x), y3 = to_fixed(p3->y);
	// check triangle winding order
	int det = det2x2(x2 - x1, x3 - x1, y2 - y1, y3 - y1);
	if (det == 0)
		return false; // zero area
	//if (det > 0)
	//	return false; // wrong winding

	// bounding box / clipping
	hs->minx = max2(fixed_ceil(min3(x1, x2, x3)), 0);
	hs->miny = max2(fixed_ceil(min3(y1, y2, y3)), 0);
	hs->maxx = min2(fixed_ceil(max3(x1, x2, x3)), SCREEN_X);
	hs->maxy = min2(fixed_ceil(max3(y1, y2, y3)), SCREEN_Y);
	if (hs->minx >= hs->maxx || hs->miny >= hs->maxy)
		return false;

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
	return true;
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


// --------------------------------------------------------------------------
// Dither3D based on https://github.com/runevision/Dither3D/blob/main/Assets/Dither3D/Dither3DInclude.cginc
// We use 4x4 dither pattern, but reduced texture XY resolution
// And a bunch of other simplifications


// equivalent to `x / exp2f((float)i)`, provided we are not in
// infinities / subnormals territory.
static inline float adjust_float_exp(float x, int i)
{
	union {
		float f;
		uint32_t u;
	} fu;
	fu.f = x;
	fu.u -= (uint32_t)i << 23;
	return fu.f;
}

// Note: no RADIAL_COMPENSATION
#define DITHER_RES (32) // Note: upstream used 64
#define DITHER_RES_SHIFT (5)
#define DITHER_RES_MASK (DITHER_RES-1)
#define DITHER_DOTS_PER_SIDE (4) // Note: upstream used (DITHER_RES/16)
#define DITHER_SLICES (DITHER_DOTS_PER_SIDE * DITHER_DOTS_PER_SIDE)

#define DITHER_SCALE (6.0f)
#define DITHER_SCALE_EXP (64.0f) // exp2f(DITHER_SCALE)

FORCE_INLINE bool do_pixel_sample(const float u, const float v, const int patternScaleLevel_i, const int subLayer_offset, const int compare, const int x, const int y, const float debug_val)
{
	// Get the UV coordinates in the current fractal level.
	// const float scaleLevelMul = 1.0f / exp2f(patternScaleLevel);
	// float uu = u * scaleLevelMul;
	// float vv = v * scaleLevelMul;
	// instead of above, we can directly alter float exponent bits:
	const float uu = adjust_float_exp(u, patternScaleLevel_i);
	const float vv = adjust_float_exp(v, patternScaleLevel_i);
	// Sample the 3D texture.
	const int u3d = (unsigned)((int)(uu * DITHER_RES)) & DITHER_RES_MASK;
	const int v3d = (unsigned)((int)(vv * DITHER_RES)) & DITHER_RES_MASK;
	const int texel_idx = (subLayer_offset + v3d) * DITHER_RES + u3d;
	const uint8_t pattern = s_dither4x4_r[texel_idx];

#if 0
	uint8_t* debug_output = plat_gfx_get_debug_frame();
	if (debug_output) {
		float3 dbg = (float3){ 0,0,0 };

		// UV
		dbg.x = fract(u), dbg.y = fract(v);

		// fractal UV
		dbg.x = fract(uu), dbg.y = fract(vv);

		// pattern
		//dbg.x = dbg.y = dbg.z = pattern / 255.0f;

		// bw
		//dbg.x = dbg.y = dbg.z = pattern < compare ? 0.0f : 1.0f;

		//dbg.x = debug_val;
		//dbg.y = -debug_val;
		//dbg.z = 0;

		const int offset = (y * SCREEN_X + x) * 4 + 0;
		debug_output[offset + 0] = (uint8_t)(saturate(dbg.x) * 255.0f);
		debug_output[offset + 1] = (uint8_t)(saturate(dbg.y) * 255.0f);
		debug_output[offset + 2] = (uint8_t)(saturate(dbg.z) * 255.0f);
		debug_output[offset + 3] = 255;
	}
#endif

	// Get the pattern value relative to the threshold, scale it
	// according to the contrast, and add the base value.
	// Note: the terms are slightly reassociated compared to upstream, since
	// we only need to threshold the result.
	return pattern < compare;
}

FORCE_INLINE bool do_pixel_sample_16_16(fixed16_16 u, fixed16_16 v, const int patternScaleLevel_i, const int subLayer_offset, const int compare, const int x, const int y, const float debug_val)
{
	// Get the UV coordinates in the current fractal level.
	if (patternScaleLevel_i >= 0) {
		u >>= patternScaleLevel_i;
		v >>= patternScaleLevel_i;
	}
	else {
		u <<= -patternScaleLevel_i;
		v <<= -patternScaleLevel_i;
	}

	// Sample the 3D texture
	const int u3d = (u >> (16 - DITHER_RES_SHIFT)) & DITHER_RES_MASK;
	const int v3d = (v >> (16 - DITHER_RES_SHIFT)) & DITHER_RES_MASK;
	const int texel_idx = (subLayer_offset + v3d) * DITHER_RES + u3d;
	const uint8_t pattern = s_dither4x4_r[texel_idx];

#if 0
	uint8_t* debug_output = plat_gfx_get_debug_frame();
	if (debug_output) {
		float3 dbg = (float3){ 0,0,0 };

		// UV
		//dbg.x = fract(u), dbg.y = fract(v);

		// fractal UV
		//dbg.x = fract(uu), dbg.y = fract(vv);
		dbg.x = u3d / (float)DITHER_RES;
		dbg.y = v3d / (float)DITHER_RES;

		// pattern
		//dbg.x = dbg.y = dbg.z = pattern / 255.0f;

		// bw
		//dbg.x = dbg.y = dbg.z = pattern < compare ? 0.0f : 1.0f;

		//dbg.x = debug_val;
		//dbg.y = -debug_val;
		//dbg.z = 0;

		const int offset = (y * SCREEN_X + x) * 4 + 0;
		debug_output[offset + 0] = (uint8_t)(saturate(dbg.x) * 255.0f);
		debug_output[offset + 1] = (uint8_t)(saturate(dbg.y) * 255.0f);
		debug_output[offset + 2] = (uint8_t)(saturate(dbg.z) * 255.0f);
		debug_output[offset + 3] = 255;
	}
#endif

	// Get the pattern value relative to the threshold, scale it
	// according to the contrast, and add the base value.
	// Note: the terms are slightly reassociated compared to upstream, since
	// we only need to threshold the result.
	return pattern < compare;
}

FORCE_INLINE int get_dither3d_compare_val(const uint8_t brightness, const uint8_t brightnessCurve, const float brightnessSpacingMultiplier)
{
	// We create sharp dots from them by increasing the contrast.
	const float _Contrast = 1.0f;
	const float contrast = _Contrast * DITHER_SCALE_EXP * brightnessSpacingMultiplier * 0.1f;
	// The base brightness value that we scale the contrast around
	// should normally be 0.5, but if the pattern is very blurred,
	// that would just make the brightness everywhere close to 0.5.
	float baseVal = lerp(128.0f, brightness, saturate(1.05f / (1.0f + contrast))) - 128.0f;
	// The brighter output we want, the lower threshold we need to use
	const uint8_t threshold = 255 - brightnessCurve;
	const int compare_val = (int)(threshold - baseVal / contrast + 0.5f);
	return compare_val;
}

FORCE_INLINE int get_dither3d_level_fraction(float spacing, int *r_patternScaleLevel_i)
{
	// Find the power-of-two level that corresponds to the dot spacing.
	//float spacingLog = log2f(spacing);
	//const float patternScaleLevel = floorf(spacingLog); // Fractal level.
	//const int patternScaleLevel_i = (int)patternScaleLevel;
	//float f = spacingLog - patternScaleLevel; // Fractional part.
	//
	// instead of above, work on float bits directly:
	union {
		float f;
		uint32_t u;
	} fu;
	fu.f = spacing;
	// patternScaleLevel is just float exponent:
	*r_patternScaleLevel_i = (int)((fu.u >> 23) & 0xFF) - 127;

	// fractional part is:
	// - take the mantissa bits of spacing,
	// - set exponent to 127, i.e. range [0,1)
	// - use that as a float and subtract 1.0
	fu.u = (fu.u & 0x7FFFFF) | 0x3F800000;
	const float f = fu.f - 1.0f;

	// Get the third coordinate for the 3D texture lookup.
	// Note: simplified the math terms compared to upstream.
	const int subLayer_i = (int)(0.75f * DITHER_SLICES * f + 0.5f);
	return subLayer_i * DITHER_RES;
}

FORCE_INLINE int get_dither3d_level_fraction_16_16(fixed16_16 spacing, int* r_patternScaleLevel_i)
{
	// As above, but in fixed point:
	// - the scale level (log2 of input) is index of highest set bit, positive or negative compared to 16
	// - fraction is the bits after that

	int level = highest_bit(spacing) - 16;
	fixed16_16 fraction = (level >= 0) ? (spacing >> level) : (spacing << -level);
	fraction &= 0xFFFF;
	*r_patternScaleLevel_i = level;

	fraction = fraction * 3 / 4 + 0x800;
	int subLayer = fraction >> 12;
	return subLayer * DITHER_RES;
}


void draw_triangle_dither3d_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6], const uint8_t brightness)
{
	raster_halfspace_t hs;
	if (!raster_halfspace_begin(&hs, p1, p2, p3, uvs))
		return;

	// Lookup brightness to make dither output have correct output
	// brightness at different input brightness values.
	const uint8_t brightnessCurve = s_dither4x4_g[brightness / 4];
	// Note: _SizeVariability fixed to default 0.0, following math simplified
	const float brightnessSpacingMultiplier = 255.0f / (brightnessCurve * 2.0f + 0.001f);
	// Scale the spacing by the specified input (power of two) scale.
	// We keep the spacing the same regardless of whether we're using
	// a pattern with more or less dots in it.
	const float spacingMul = DITHER_SCALE_EXP * DITHER_DOTS_PER_SIDE * 0.125f * brightnessSpacingMultiplier;
	const int compare_val = get_dither3d_compare_val(brightness, brightnessCurve, brightnessSpacingMultiplier);

	// Rasterize: rows
	raster_halfspace_state_t state;
	for (raster_halfspace_y_begin(&hs, &state); raster_halfspace_y_continue(&hs, &state); raster_halfspace_y_step(&hs, &state))
	{
		uint8_t* output_0 = bitmap + state.y * rowstride;
		uint8_t* output_1 = output_0 + rowstride;		

		// Rasterize: columns
		for (raster_halfspace_x_begin(&hs, &state); raster_halfspace_x_continue(&hs, &state); raster_halfspace_x_step(&hs, &state))
		{
			int inner = raster_halfspace_x_inner(&hs, &state);
			if (inner == 1) break;
			if (inner == 2) continue;

			// We define a spacing variable which linearly correlates with
			// the average distance between dots.
			// Note: simpler than upstream which calculates two frequency values
			// based on singular value decomposition of derivatives.
			float spacing = (fabsf(state.dudx) + fabsf(state.dvdx) + fabsf(state.dudy) + fabsf(state.dvdy)) * 0.25f;
			spacing *= spacingMul;

			int patternScaleLevel_i;
			const int subLayer_offset = get_dither3d_level_fraction(spacing, &patternScaleLevel_i);

			// Note: accumulating pixel outputs/masks into 32 bit words and writing them out
			// once they are done (or trailing after the loop) seems to be slightly slower on Playdate,
			// than simply operating on bytes in a naive way.
			const int byte_idx = state.x / 8;
			const int mask_0 = 1 << (7 - ((state.x + 0) & 7));
			const int mask_1 = 1 << (7 - ((state.x + 1) & 7));

			float dbg = 0.0f;
			//dbg = state.dudx * 15;
			//dbg = state.dvdx * 15;
			//dbg = state.dudy * 15;
			//dbg = state.dvdy * 15;
			//dbg = spacing;
			dbg = patternScaleLevel_i / 16.0f;
			dbg = subLayer_offset / (float)(DITHER_RES * DITHER_SLICES);
			if (state.inside_00)
			{
				bool check = do_pixel_sample(state.u_00, state.v_00, patternScaleLevel_i, subLayer_offset, compare_val, state.x, state.y, dbg);
				if (check)
					output_0[byte_idx] &= ~mask_0;
				else
					output_0[byte_idx] |= mask_0;
			}
			if (state.inside_01)
			{
				bool check = do_pixel_sample(state.u_01, state.v_01, patternScaleLevel_i, subLayer_offset, compare_val, state.x + 1, state.y, dbg);
				if (check)
					output_0[byte_idx] &= ~mask_1;
				else
					output_0[byte_idx] |= mask_1;
			}
			if (state.inside_10)
			{
				bool check = do_pixel_sample(state.u_10, state.v_10, patternScaleLevel_i, subLayer_offset, compare_val, state.x, state.y + 1, dbg);
				if (check)
					output_1[byte_idx] &= ~mask_0;
				else
					output_1[byte_idx] |= mask_0;
			}
			if (state.inside_11)
			{
				bool check = do_pixel_sample(state.u_11, state.v_11, patternScaleLevel_i, subLayer_offset, compare_val, state.x + 1, state.y + 1, dbg);
				if (check)
					output_1[byte_idx] &= ~mask_1;
				else
					output_1[byte_idx] |= mask_1;
			}
		}
	}
}

typedef struct spacing_gradients {
	float sz[3]; // spacing/z for each vertex
	float sz_dx, sz_dy;
} spacing_gradients;

typedef struct spacing_edge {
	float sz, sz_step, sz_step_extra;
} spacing_edge;

static void draw_scanline_dither3d(uint8_t* bitmap, int rowstride,
	const tri_gradients* grad, const tri_edge_uv* edge_l, const tri_edge_uv* edge_r,
	const spacing_edge* spc_l, const spacing_edge* spc_r,
	const uint8_t compare_val, const float spacingMul)
{
	raster_scanline_line_t line;
	if (!raster_scanline_line_init(&line, grad, edge_l, edge_r))
		return;

	int X = line.x_start;
	bitmap += edge_l->Y * rowstride;
	// write out pixels 32 at a time
	uint32_t mask = 0;
	uint32_t* p = (uint32_t*)bitmap + X / 32;
	uint32_t color = 0;


	const fixed16_16 spacing_l = FloatToFixed16_16(spc_l->sz / edge_l->invz);
	const fixed16_16 spacing_r = FloatToFixed16_16(spc_r->sz / edge_r->invz);
	const fixed16_16 spacing_step = (spacing_r - spacing_l) / (line.x_end - line.x_start);
	fixed16_16 spacing = spacing_l;

	while (raster_scanline_spans_continue(&line))
	{
		raster_scanline_inner_begin(&line);

		for (int i = 0; i < kAffineLength; i++)
		{
			int patternScaleLevel_i;
			const int subLayer_offset = get_dither3d_level_fraction_16_16(spacing, &patternScaleLevel_i);
			float dbg = 0.0f;
			//dbg = Fixed16_16ToFloat(spacing);
			//dbg = patternScaleLevel_i / 16.0f;
			//dbg = subLayer_offset / (float)(DITHER_RES * DITHER_SLICES);
			bool check = do_pixel_sample_16_16(line.U, line.V, patternScaleLevel_i, subLayer_offset, compare_val, X, edge_l->Y, dbg);

			//int bit_mask = 1 << (7 - (X & 7));
			//if (check)
			//	bitmap[X / 8] &= ~bit_mask;
			//else
			//	bitmap[X / 8] |= bit_mask;
			mask |= 0x80000000u >> (X & 31);
			uint32_t texpix = check ? 0x0 : 0x80;
			color |= (texpix << 24) >> (X % 32);


			X++;
			spacing += spacing_step;
			raster_scanline_inner_step(&line);

			if (X % 32 == 0)
			{
				_drawMaskPattern(p++, swap(mask), swap(color));
				mask = 0;
				color = 0;
			}
		}
		raster_scanline_spans_step(&line);
	}
	if (raster_scanline_has_rem(&line))
	{
		raster_scanline_rem_begin(&line, grad, edge_r);
		for (int i = 0; i <= line.aff_remainder; i++)
		{
			int patternScaleLevel_i;
			const int subLayer_offset = get_dither3d_level_fraction_16_16(spacing, &patternScaleLevel_i);
			float dbg = 0.0f;
			//dbg = Fixed16_16ToFloat(spacing);
			//dbg = patternScaleLevel_i / 16.0f;
			//dbg = subLayer_offset / (float)(DITHER_RES * DITHER_SLICES);
			bool check = do_pixel_sample_16_16(line.U, line.V, patternScaleLevel_i, subLayer_offset, compare_val, X, edge_l->Y, dbg);

			//int bit_mask = 1 << (7 - (X & 7));
			//if (check)
			//	bitmap[X / 8] &= ~bit_mask;
			//else
			//	bitmap[X / 8] |= bit_mask;
			mask |= 0x80000000u >> (X & 31);
			uint32_t texpix = check ? 0x0 : 0x80;
			color |= (texpix << 24) >> (X % 32);

			X++;
			spacing += spacing_step;
			raster_scanline_inner_step(&line);

			if (X % 32 == 0)
			{
				_drawMaskPattern(p++, swap(mask), swap(color));
				mask = 0;
				color = 0;
			}
		}
	}

	_drawMaskPattern(p, swap(mask), swap(color));
}

static void spacing_gradients_init(spacing_gradients* t, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const tri_gradients* grad, const float spacingMul)
{
	for (int i = 0; i < 3; ++i)
	{
		const float ux = (grad->uz[i] + grad->uz_dx) / (grad->invz[i] + grad->invz_dx);
		const float vx = (grad->vz[i] + grad->vz_dx) / (grad->invz[i] + grad->invz_dx);
		const float uy = (grad->uz[i] + grad->uz_dy) / (grad->invz[i] + grad->invz_dy);
		const float vy = (grad->vz[i] + grad->vz_dy) / (grad->invz[i] + grad->invz_dy);
		const float dudx = ux - uvs[i * 2 + 0];
		const float dvdx = vx - uvs[i * 2 + 1];
		const float dudy = uy - uvs[i * 2 + 0];
		const float dvdy = vy - uvs[i * 2 + 1];
		float spacing = (fabsf(dudx) + fabsf(dvdx) + fabsf(dudy) + fabsf(dvdy)) * 0.25f;
		spacing *= spacingMul;
		t->sz[i] = spacing * grad->invz[i];
	}

	float invdx = grad->inv_det;
	float invdy = -grad->inv_det;

	float x02 = p0->x - p2->x;
	float y02 = p0->y - p2->y;
	float x12 = p1->x - p2->x;
	float y12 = p1->y - p2->y;
	t->sz_dx = invdx * (((t->sz[1] - t->sz[2]) * y02) - ((t->sz[0] - t->sz[2]) * y12));
	t->sz_dy = invdy * (((t->sz[1] - t->sz[2]) * x02) - ((t->sz[0] - t->sz[2]) * x12));
}

static void spacing_edge_init(spacing_edge* t, const tri_edge_uv* edge, const spacing_gradients* grad, int top)
{
	t->sz = grad->sz[top] + edge->y_prestep * grad->sz_dy + edge->x_prestep * grad->sz_dx;
	t->sz_step = edge->x_step * grad->sz_dx + grad->sz_dy;
	t->sz_step_extra = grad->sz_dx;
}

FORCE_INLINE void spacing_edge_step(tri_edge_uv* t, spacing_edge* s)
{
	t->X += t->x_step; t->Y++; t->height--;
	t->uz += t->uz_step; t->vz += t->vz_step; t->invz += t->invz_step;
	s->sz += s->sz_step;

	t->error_term += t->numerator;
	if (t->error_term >= t->denominator) {
		t->X++;
		t->error_term -= t->denominator;
		t->invz += t->invz_step_extra;
		t->uz += t->uz_step_extra; t->vz += t->vz_step_extra;
		s->sz += s->sz_step_extra;
	}
}

static void draw_triangle_dither3d_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const uint8_t brightness)
{
	raster_scanline_uv_t tri;
	if (!raster_scanline_uv_begin(&tri, p0, p1, p2, uvs, 1.0f))
		return;

	// Lookup brightness to make dither output have correct output
	// brightness at different input brightness values.
	const uint8_t brightnessCurve = s_dither4x4_g[brightness / 4];
	// Note: _SizeVariability fixed to default 0.0, following math simplified
	const float brightnessSpacingMultiplier = 255.0f / (brightnessCurve * 2.0f + 0.001f);
	// Scale the spacing by the specified input (power of two) scale.
	// We keep the spacing the same regardless of whether we're using
	// a pattern with more or less dots in it.
	const float spacingMul = DITHER_SCALE_EXP * DITHER_DOTS_PER_SIDE * 0.125f * brightnessSpacingMultiplier;
	const int compare_val = get_dither3d_compare_val(brightness, brightnessCurve, brightnessSpacingMultiplier);

	// Calculate spacing at vertices of triangle, interpolate across
	spacing_gradients s_grad;
	spacing_gradients_init(&s_grad, p0, p1, p2, uvs, &tri.grad, spacingMul);
	spacing_edge s_top_bottom, s_top_mid, s_mid_bot;
	spacing_edge_init(&s_top_bottom, &tri.e_top_bottom, &s_grad, tri.v_top);
	spacing_edge_init(&s_top_mid, &tri.e_top_mid, &s_grad, tri.v_top);
	spacing_edge_init(&s_mid_bot, &tri.e_mid_bottom, &s_grad, tri.v_mid);
	spacing_edge* s_e_left;
	spacing_edge* s_e_right;

	int height = raster_scanline_uv_begin1(&tri);
	s_e_left = tri.mid_is_left ? &s_top_bottom : &s_top_mid;
	s_e_right = tri.mid_is_left ? &s_top_mid : &s_top_bottom;
	while (height--) {
		draw_scanline_dither3d(bitmap, rowstride, &tri.grad, tri.e_left, tri.e_right, s_e_left, s_e_right, compare_val, spacingMul);
		spacing_edge_step(tri.e_left, s_e_left);
		spacing_edge_step(tri.e_right, s_e_right);
	}
	height = raster_scanline_uv_begin2(&tri);
	s_e_left = tri.mid_is_left ? &s_top_bottom : &s_mid_bot;
	s_e_right = tri.mid_is_left ? &s_mid_bot : &s_top_bottom;
	while (height--) {
		draw_scanline_dither3d(bitmap, rowstride, &tri.grad, tri.e_left, tri.e_right, s_e_left, s_e_right, compare_val, spacingMul);
		spacing_edge_step(tri.e_left, s_e_left);
		spacing_edge_step(tri.e_right, s_e_right);
	}
}

// --------------------------------------------------------------------------
// Simple pure black/white checkerboard

#define DEBUG_CHECKER_RENDER 0

static void draw_scanline_checker(uint8_t* bitmap, int rowstride, const tri_gradients* grad, const tri_edge_uv* edge_l, const tri_edge_uv* edge_r)
{
	raster_scanline_line_t line;
	if (!raster_scanline_line_init(&line, grad, edge_l, edge_r))
		return;

	int X = line.x_start;

	bitmap += edge_l->Y * rowstride;
	// write out pixels 32 at a time
	uint32_t mask = 0;
	uint32_t* p = (uint32_t*)bitmap + X / 32;
	uint32_t color = 0;

	while (raster_scanline_spans_continue(&line))
	{
		raster_scanline_inner_begin(&line);
		for (int i = 0; i < kAffineLength; i++)
		{
			bool checker = ((line.U & 0xFFFF) >= 0x8000) != ((line.V & 0xFFFF) >= 0x8000);
			//checker = true;

			mask |= 0x80000000u >> (X & 31);
			uint32_t texpix = checker ? 0x0 : 0x80;
			color |= (texpix << 24) >> (X % 32);
			//int bit_mask = 1 << (7 - (X & 7));
			//if (checker)
			//	bitmap[X / 8] &= ~bit_mask;
			//else
			//	bitmap[X / 8] |= bit_mask;

#if DEBUG_CHECKER_RENDER
			uint8_t* dbg = plat_gfx_get_debug_frame();
			if (dbg)
			{
				int dbg_idx = (edge_l->Y * SCREEN_X + X) * 4;
				dbg[dbg_idx + 0] += checker ? 250 : 0;
				dbg[dbg_idx + 1] += (line.U / 5 >> 8) & 0xFF;
				dbg[dbg_idx + 2] += (line.V / 5 >> 8) & 0xFF;
				dbg[dbg_idx + 3] = 255;
			}
#endif

			X++;
			raster_scanline_inner_step(&line);

			if (X % 32 == 0)
			{
				_drawMaskPattern(p++, swap(mask), swap(color));
				mask = 0;
				color = 0;
			}
		}
		raster_scanline_spans_step(&line);
	}

	if (raster_scanline_has_rem(&line))
	{
		raster_scanline_rem_begin(&line, grad, edge_r);
		for (int i = 0; i <= line.aff_remainder; i++)
		{
			bool checker = ((line.U & 0xFFFF) >= 0x8000) != ((line.V & 0xFFFF) >= 0x8000);
			//checker = true;

			mask |= 0x80000000u >> (X & 31);
			uint32_t texpix = checker ? 0x0 : 0x80;
			color |= (texpix << 24) >> (X % 32);

			//int bit_mask = 1 << (7 - (X & 7));
			//if (checker)
			//	bitmap[X / 8] &= ~bit_mask;
			//else
			//	bitmap[X / 8] |= bit_mask;
#if DEBUG_CHECKER_RENDER
			uint8_t* dbg = plat_gfx_get_debug_frame();
			if (dbg)
			{
				int dbg_idx = (edge_l->Y * SCREEN_X + X) * 4;
				dbg[dbg_idx + 0] += checker ? 250 : 0;
				dbg[dbg_idx + 1] += (line.U / 5 >> 8) & 0xFF;
				dbg[dbg_idx + 2] += (line.V / 5 >> 8) & 0xFF;
				dbg[dbg_idx + 3] = 255;
			}
#endif

			X++;
			raster_scanline_inner_step(&line);

			if (X % 32 == 0)
			{
				_drawMaskPattern(p++, swap(mask), swap(color));
				mask = 0;
				color = 0;
			}
		}
	}

	_drawMaskPattern(p, swap(mask), swap(color));
}

static void draw_triangle_checker_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6])
{
	raster_scanline_uv_t tri;
	if (!raster_scanline_uv_begin(&tri, p0, p1, p2, uvs, 5.0f))
		return;

	int height = raster_scanline_uv_begin1(&tri);
	while (height--) {
		draw_scanline_checker(bitmap, rowstride, &tri.grad, tri.e_left, tri.e_right);
		raster_scanline_uv_y_step(&tri);
	}
	height = raster_scanline_uv_begin2(&tri);
	while (height--) {
		draw_scanline_checker(bitmap, rowstride, &tri.grad, tri.e_left, tri.e_right);
		raster_scanline_uv_y_step(&tri);
	}
}

FORCE_INLINE bool do_checker(float u, float v, int x, int y)
{
	float uu = fract(u * 5.0f);
	float vv = fract(v * 5.0f);
	bool checker = (uu > 0.5f) != (vv > 0.5f);
	//checker = true;

#if DEBUG_CHECKER_RENDER
	uint8_t* dbg = plat_gfx_get_debug_frame();
	if (dbg)
	{
		int dbg_idx = (y * SCREEN_X + x) * 4;

		int ui = ((int)(fract(u) * 64)) & 63;
		int vi = ((int)(fract(v) * 64)) & 63;
		dbg[dbg_idx + 0] += checker ? 250 : 0;
		dbg[dbg_idx + 1] += ui * 4;
		dbg[dbg_idx + 2] += vi * 4;

		dbg[dbg_idx + 3] = 255;
	}
#endif

	return checker;
}

void draw_triangle_checker_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6])
{
	raster_halfspace_t hs;
	if (!raster_halfspace_begin(&hs, p1, p2, p3, uvs))
		return;

	raster_halfspace_state_t state;
	for (raster_halfspace_y_begin(&hs, &state); raster_halfspace_y_continue(&hs, &state); raster_halfspace_y_step(&hs, &state))
	{
		uint8_t* output_0 = bitmap + state.y * rowstride;
		uint8_t* output_1 = output_0 + rowstride;
		for (raster_halfspace_x_begin(&hs, &state); raster_halfspace_x_continue(&hs, &state); raster_halfspace_x_step(&hs, &state))
		{
			int inner = raster_halfspace_x_inner(&hs, &state);
			if (inner == 1) break;
			if (inner == 2) continue;

			const int byte_idx = state.x / 8;
			const int mask_0 = 1 << (7 - ((state.x + 0) & 7));
			const int mask_1 = 1 << (7 - ((state.x + 1) & 7));

			if (state.inside_00)
			{
				bool check = do_checker(state.u_00, state.v_00, state.x, state.y);
				if (check)
					output_0[byte_idx] &= ~mask_0;
				else
					output_0[byte_idx] |= mask_0;
			}
			if (state.inside_01)
			{
				bool check = do_checker(state.u_01, state.v_01, state.x + 1, state.y);
				if (check)
					output_0[byte_idx] &= ~mask_1;
				else
					output_0[byte_idx] |= mask_1;
			}
			if (state.inside_10)
			{
				bool check = do_checker(state.u_10, state.v_10, state.x, state.y + 1);
				if (check)
					output_1[byte_idx] &= ~mask_0;
				else
					output_1[byte_idx] |= mask_0;
			}
			if (state.inside_11)
			{
				bool check = do_checker(state.u_11, state.v_11, state.x + 1, state.y + 1);
				if (check)
					output_1[byte_idx] &= ~mask_1;
				else
					output_1[byte_idx] |= mask_1;
			}
		}
	}
}

// --------------------------------------------------------------------------
// Overall scene rendering


void scene_init(Scene* scene)
{
	scene_setCamera(scene, (float3) { 0, 0, 0 }, (float3) { 0, 0, 1 }, 1.0, (float3) { 0, 1, 0 });
	scene_setLight(scene, (float3) { 0, -1, 0 });

	scene_setCenter(scene, 0.5, 0.5);

	scene->tmp_points_cap = scene->tmp_faces_cap = 0;
	scene->tmp_points = NULL;
	scene->tmp_face_normals = NULL;
	scene->tmp_order_table = NULL;
	scene->tmp_dist_table = NULL;
}

void scene_shutdown(Scene* scene)
{
	plat_free(scene->tmp_points);
	plat_free(scene->tmp_face_normals);
	plat_free(scene->tmp_order_table);
	plat_free(scene->tmp_dist_table);
}

void scene_setLight(Scene* scene, float3 light)
{
	scene->light = v3_normalize(light);
}

void scene_setCenter(Scene* scene, float x, float y)
{
	scene->centerx = x;
	scene->centery = y;
}

void scene_setCamera(Scene* scene, float3 origin, float3 lookAt, float scale, float3 up)
{
	xform camera = xform_identity;

	camera.x = -origin.x;
	camera.y = -origin.y;
	camera.z = -origin.z;

	float3 dir = (float3){ lookAt.x - origin.x, lookAt.y - origin.y, lookAt.z - origin.z };

	float l = sqrtf(v3_lensq(&dir));

	dir.x /= l;
	dir.y /= l;
	dir.z /= l;

	scene->scale = 240 * scale;

	// first yaw around the y axis

	float h = 0;

	if (dir.x != 0 || dir.z != 0)
	{
		h = sqrtf(dir.x * dir.x + dir.z * dir.z);

		xform yaw = xform_make(dir.z / h, 0, -dir.x / h, 0, 1, 0, dir.x / h, 0, dir.z / h);
		camera = xform_multiply(&camera, &yaw);
	}

	// then pitch up/down to y elevation

	xform pitch = xform_make(1, 0, 0, 0, h, -dir.y, 0, dir.y, h);
	camera = xform_multiply(&camera, &pitch);

	// and roll to position the up vector

	if (up.x != 0 || up.y != 0)
	{
		l = sqrtf(up.x * up.x + up.y * up.y);
		xform roll = xform_make(up.y / l, up.x / l, 0, -up.x / l, up.y / l, 0, 0, 0, 1);

		scene->camera = xform_multiply(&camera, &roll);
	}
	else
		scene->camera = camera;
}

typedef uint8_t Pattern[8];

static Pattern patterns[] =
{
	{ 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 },
	{ 0x80, 0x00, 0x00, 0x00, 0x08, 0x00, 0x00, 0x00 },
	{ 0x88, 0x00, 0x00, 0x00, 0x88, 0x00, 0x00, 0x00 },
	{ 0x88, 0x00, 0x20, 0x00, 0x88, 0x00, 0x02, 0x00 },
	{ 0x88, 0x00, 0x22, 0x00, 0x88, 0x00, 0x22, 0x00 },
	{ 0xa8, 0x00, 0x22, 0x00, 0x8a, 0x00, 0x22, 0x00 },
	{ 0xaa, 0x00, 0x22, 0x00, 0xaa, 0x00, 0x22, 0x00 },
	{ 0xaa, 0x00, 0xa2, 0x00, 0xaa, 0x00, 0x2a, 0x00 },
	{ 0xaa, 0x00, 0xaa, 0x00, 0xaa, 0x00, 0xaa, 0x00 },
	{ 0xaa, 0x40, 0xaa, 0x00, 0xaa, 0x04, 0xaa, 0x00 },
	{ 0xaa, 0x44, 0xaa, 0x00, 0xaa, 0x44, 0xaa, 0x00 },
	{ 0xaa, 0x44, 0xaa, 0x10, 0xaa, 0x44, 0xaa, 0x01 },
	{ 0xaa, 0x44, 0xaa, 0x11, 0xaa, 0x44, 0xaa, 0x11 },
	{ 0xaa, 0x54, 0xaa, 0x11, 0xaa, 0x45, 0xaa, 0x11 },
	{ 0xaa, 0x55, 0xaa, 0x11, 0xaa, 0x55, 0xaa, 0x11 },
	{ 0xaa, 0x55, 0xaa, 0x51, 0xaa, 0x55, 0xaa, 0x15 },
	{ 0xaa, 0x55, 0xaa, 0x55, 0xaa, 0x55, 0xaa, 0x55 },
	{ 0xba, 0x55, 0xaa, 0x55, 0xab, 0x55, 0xaa, 0x55 },
	{ 0xbb, 0x55, 0xaa, 0x55, 0xbb, 0x55, 0xaa, 0x55 },
	{ 0xbb, 0x55, 0xea, 0x55, 0xbb, 0x55, 0xae, 0x55 },
	{ 0xbb, 0x55, 0xee, 0x55, 0xbb, 0x55, 0xee, 0x55 },
	{ 0xfb, 0x55, 0xee, 0x55, 0xbf, 0x55, 0xee, 0x55 },
	{ 0xff, 0x55, 0xee, 0x55, 0xff, 0x55, 0xee, 0x55 },
	{ 0xff, 0x55, 0xfe, 0x55, 0xff, 0x55, 0xef, 0x55 },
	{ 0xff, 0x55, 0xff, 0x55, 0xff, 0x55, 0xff, 0x55 },
	{ 0xff, 0x55, 0xff, 0xd5, 0xff, 0x55, 0xff, 0x5d },
	{ 0xff, 0x55, 0xff, 0xdd, 0xff, 0x55, 0xff, 0xdd },
	{ 0xff, 0x75, 0xff, 0xdd, 0xff, 0x57, 0xff, 0xdd },
	{ 0xff, 0x77, 0xff, 0xdd, 0xff, 0x77, 0xff, 0xdd },
	{ 0xff, 0x77, 0xff, 0xfd, 0xff, 0x77, 0xff, 0xdf },
	{ 0xff, 0x77, 0xff, 0xff, 0xff, 0x77, 0xff, 0xff },
	{ 0xff, 0xf7, 0xff, 0xff, 0xff, 0x7f, 0xff, 0xff },
	{ 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff }
};

void drawShapeFace(const Scene* scene, uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float3* normal, const Mesh* mesh, int tri_index, enum DrawStyle style, bool wire)
{
	// lighting
	float v = v3_dot(*normal, scene->light) * 0.5f + 0.5f;
	//v = normal->z * 0.5f + 0.5f;

	// cheap gamma adjust
	//v = v * v;

	if (style == Draw_Pattern)
	{
		int vi = (int)(32.99f * v);
		if (vi > 32)
			vi = 32;
		else if (vi < 0)
			vi = 0;

		const uint8_t* pattern = (const uint8_t*)&patterns[vi];
		draw_triangle_pattern_scanline(bitmap, rowstride, p1, p2, p3, pattern);
	}
	else if (style == Draw_Bluenoise)
	{
		// draw into byte buffer for blue noise thresholding
		int col = (int)(v * 255.0f);
		if (col < 0) col = 0;
		if (col > 255) col = 255;

		draw_triangle_bluenoise_scanline(bitmap, rowstride, p1, p2, p3, col);
	}
	else if (style == Draw_Checker_Scanline)
	{
		// draw strictly black/white checker based on UV coordinates
		draw_triangle_checker_scanline(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6);
	}
	else if (style == Draw_Checker_Halfspace)
	{
		// draw strictly black/white checker based on UV coordinates
		draw_triangle_checker_halfspace(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6);
	}
	else if (style == Draw_Dither3D_Scanline)
	{
		draw_triangle_dither3d_scanline(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6, (uint8_t)(saturate(v) * 255.0f));
	}
	else if (style == Draw_Dither3D_Halfspace)
	{
		// draw strictly black/white checker based on UV coordinates
		draw_triangle_dither3d_halfspace(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6, (uint8_t)(saturate(v) * 255.0f));
	}

	// wireframe
	if (wire)
	{
		const uint8_t* pattern = (const uint8_t*)&patterns[0];
		drawLine(bitmap, rowstride, p1, p2, 1, pattern);
		drawLine(bitmap, rowstride, p2, p3, 1, pattern);
		drawLine(bitmap, rowstride, p3, p1, 1, pattern);
	}
}


static Scene* s_facesort_instance;
static int compareFaceZ(const void* a, const void* b)
{
	uint16_t idxa = *(const uint16_t*)a;
	uint16_t idxb = *(const uint16_t*)b;
	float za = s_facesort_instance->tmp_dist_table[idxa];
	float zb = s_facesort_instance->tmp_dist_table[idxb];
	if (za < zb)
		return +1;
	if (za > zb)
		return -1;
	return 0;
}

void scene_drawMesh(Scene* scene, uint8_t* buffer, int rowstride, const Mesh* mesh, const xform* matrix, enum DrawStyle style, bool wire)
{
	// temporary buffers
	if (scene->tmp_points_cap < mesh->vertex_count) {
		scene->tmp_points_cap = mesh->vertex_count;
		scene->tmp_points = plat_realloc(scene->tmp_points, scene->tmp_points_cap * sizeof(scene->tmp_points[0]));
	}
	if (scene->tmp_faces_cap < mesh->tri_count) {
		scene->tmp_faces_cap = mesh->tri_count;
		scene->tmp_face_normals = plat_realloc(scene->tmp_face_normals, scene->tmp_faces_cap * sizeof(scene->tmp_face_normals[0]));
		scene->tmp_order_table = plat_realloc(scene->tmp_order_table, scene->tmp_faces_cap * sizeof(scene->tmp_order_table[0]));
		scene->tmp_dist_table = plat_realloc(scene->tmp_dist_table, scene->tmp_faces_cap * sizeof(scene->tmp_dist_table[0]));
	}

	// transform points
	for (int i = 0; i < mesh->vertex_count; ++i)
		scene->tmp_points[i] = xform_transform_pt(&scene->camera, xform_transform_pt(matrix, mesh->vertices[i]));

	// compute face normals and midpoints
	const uint16_t* ibPtr = mesh->tris;
	for (int i = 0; i < mesh->tri_count; ++i, ibPtr += 3)
	{
		uint16_t idx0 = ibPtr[0];
		uint16_t idx1 = ibPtr[1];
		uint16_t idx2 = ibPtr[2];
		scene->tmp_order_table[i] = i;
		scene->tmp_face_normals[i] = v3_tri_normal(&scene->tmp_points[idx0], &scene->tmp_points[idx1], &scene->tmp_points[idx2]);

		float z = (scene->tmp_points[idx0].z + scene->tmp_points[idx1].z + scene->tmp_points[idx2].z);
		scene->tmp_dist_table[i] = z;
	}

	// project points to screen
	for (int i = 0; i < mesh->vertex_count; ++i)
	{
		float3* p = &scene->tmp_points[i];
		if (p->z > 0)
		{
			p->x = scene->scale * (p->x / p->z + 1.6666666f * scene->centerx);
			p->y = scene->scale * (p->y / p->z + scene->centery);
		}
	}

	// backface / z / screen bounds cull triangles
	int vis_tri_count = 0;
	ibPtr = mesh->tris;
	for (int i = 0; i < mesh->tri_count; ++i, ibPtr += 3)
	{
		uint16_t idx0 = ibPtr[0];
		uint16_t idx1 = ibPtr[1];
		uint16_t idx2 = ibPtr[2];
		const float3* p0 = &scene->tmp_points[idx0];
		const float3* p1 = &scene->tmp_points[idx1];
		const float3* p2 = &scene->tmp_points[idx2];

		// If any vertex is behind the camera, skip it
		if (p0->z <= 0 || p1->z <= 0 || p2->z <= 0)
			continue;

		// quick bounds check
		float x1 = p0->x;
		float y1 = p0->y;
		float x2 = p1->x;
		float y2 = p1->y;
		float x3 = p2->x;
		float y3 = p2->y;
		if ((x1 < 0 && x2 < 0 && x3 < 0) ||
			(x1 >= SCREEN_X && x2 >= SCREEN_X && x3 >= SCREEN_X) ||
			(y1 < 0 && y2 < 0 && y3 < 0) ||
			(y1 >= SCREEN_Y && y2 >= SCREEN_Y && y3 >= SCREEN_Y))
			continue;

		// only render front side of faces via winding order
		float dx21 = x2 - x1;
		float dy31 = y3 - y1;
		float dx31 = x3 - x1;
		float dy21 = y2 - y1;
		float d = dx21 * dy31 - dy21 * dx31;
		if (d >= 0)
			continue;

		scene->tmp_order_table[vis_tri_count] = i;
		++vis_tri_count;
	}

	// sort faces by z
	s_facesort_instance = scene;
	qsort(scene->tmp_order_table, vis_tri_count, sizeof(scene->tmp_order_table[0]), compareFaceZ);

	// draw faces
	for (int i = 0; i < vis_tri_count; ++i)
	{
		uint16_t fi = scene->tmp_order_table[i];
		uint16_t idx0 = mesh->tris[fi * 3 + 0];
		uint16_t idx1 = mesh->tris[fi * 3 + 1];
		uint16_t idx2 = mesh->tris[fi * 3 + 2];
		drawShapeFace(scene, buffer, rowstride, &scene->tmp_points[idx0], &scene->tmp_points[idx1], &scene->tmp_points[idx2], &scene->tmp_face_normals[fi], mesh, fi, style, wire);
	}
}
