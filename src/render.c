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
// - mini3d-plus, https://github.com/nstbayless/mini3d-plus, MIT license


static inline void
_drawMaskPattern(uint32_t* p, uint32_t mask, uint32_t color)
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
// Classical scan-line rasterizer, based on Playdate SDK mini3d-plus and Github mini3d-plus

static inline void sortTri_t(const float3** p1, const float3** p2, const float3** p3, float2* t1, float2* t2, float2* t3)
{
	float y1 = (*p1)->y, y2 = (*p2)->y, y3 = (*p3)->y;

	if (y1 <= y2 && y1 < y3)
	{
		if (y3 < y2) // 1,3,2
		{
			const float3* tmp = *p2;
			float2 tmpt = *t2;
			*p2 = *p3; *t2 = *t3;
			*p3 = tmp; *t3 = tmpt;
		}
	}
	else if (y2 < y1 && y2 < y3)
	{
		const float3* tmp = *p1;
		float2 tmpt = *t1;
		*p1 = *p2; *t1 = *t2;

		if (y3 < y1) // 2,3,1
		{
			*p2 = *p3; *t2 = *t3;
			*p3 = tmp; *t3 = tmpt;
		}
		else // 2,1,3
		{
			*p2 = tmp; *t2 = tmpt;
		}
	}
	else
	{
		const float3* tmp = *p1;
		float2 tmpt = *t1;
		*p1 = *p3; *t1 = *t3;

		if (y1 < y2) // 3,1,2
		{
			*p3 = *p2; *t3 = *t2;
			*p2 = tmp; *t2 = tmpt;
		}
		else // 3,2,1
		{
			*p3 = tmp; *t3 = tmpt;
		}
	}
}

#define UV_SHIFT 21
#define W_SHIFT 28

typedef struct raster_scanline_t
{
	// These are used in the scanline inner loop
	int x, endx; // in pixels: current position, end position of the scanline
	int u, v, w; // current U, V, W (1/Z) in fixed point
	int dudx, dvdx, dwdx; // u, v, w steps in x direction

	// These are used in the vertical loop
	int y, y2, y3; // in pixels: current position; middle/bottom rows
	bool part_12; // we are processing y1-y2 part currently
	int x1row, x2row, urow, vrow, wrow;

	int dx1, dx2, dudy, dvdy, dwdy;

	int sb, sc;

	float w1, w2, w3;
	float u1, v1, u2, v2, u3, v3;

	// Triangle vertices
	const float3* p1;
	const float3* p2;
	const float3* p3;
	// Triangle UVs
	float2 t1, t2, t3;
} raster_scanline_t;

FORCE_INLINE bool raster_scanline_begin(raster_scanline_t* hs, const float3* pos1, const float3* pos2, const float3* pos3, const float uvs[6], float texDimPowerOfTwo)
{
	hs->p1 = pos1, hs->p2 = pos2, hs->p3 = pos3;
	hs->t1 = (float2){ uvs[0], uvs[1] };
	hs->t2 = (float2){ uvs[2], uvs[3] };
	hs->t3 = (float2){ uvs[4], uvs[5] };
	// Order vertices so that p1 is at top, p2 in the middle, p3 at the bottom
	sortTri_t(&hs->p1, &hs->p2, &hs->p3, &hs->t1, &hs->t2, &hs->t3);

	const int y1 = (int)hs->p1->y;
	hs->y2 = (int)hs->p2->y;
	hs->y3 = (int)hs->p3->y;
	if (y1 >= SCREEN_Y || hs->y3 < 0)
		return false; // triangle outside of Y range

	int det = (int)((pos3->x - pos1->x) * (pos2->y - pos1->y) - (pos2->x - pos1->x) * (pos3->y - pos1->y));
	if (det == 0)
		return false; // zero area

	hs->y2 = min2(SCREEN_Y, hs->y2);
	hs->y3 = min2(SCREEN_Y, hs->y3);
	hs->y = y1;
	hs->part_12 = true;

	// scale UVs to texture size
	hs->t1.x *= texDimPowerOfTwo; hs->t1.y *= texDimPowerOfTwo;
	hs->t2.x *= texDimPowerOfTwo; hs->t2.y *= texDimPowerOfTwo;
	hs->t3.x *= texDimPowerOfTwo; hs->t3.y *= texDimPowerOfTwo;

	hs->x1row = (int)(hs->p1->x * (1 << 16));
	hs->x2row = hs->x1row;

	hs->sb = slope(hs->p1->x, hs->p1->y, hs->p2->x, hs->p2->y, 16);
	hs->sc = slope(hs->p1->x, hs->p1->y, hs->p3->x, hs->p3->y, 16);

	hs->dx1 = min2(hs->sb, hs->sc);
	hs->dx2 = max2(hs->sb, hs->sc);

	hs->w1 = 1.0f / hs->p1->z;
	hs->w2 = 1.0f / hs->p2->z;
	hs->w3 = 1.0f / hs->p3->z;
	hs->u1 = hs->t1.x * hs->w1;
	hs->v1 = hs->t1.y * hs->w1;
	hs->u2 = hs->t2.x * hs->w2;
	hs->v2 = hs->t2.y * hs->w2;
	hs->u3 = hs->t3.x * hs->w3;
	hs->v3 = hs->t3.y * hs->w3;

	const float inv_dy31 = 1.0f / (hs->p3->y - hs->p1->y);
	const float dy21 = hs->p2->y - hs->p1->y;
	const float mx = hs->p1->x + dy21 * (hs->p3->x - hs->p1->x) * inv_dy31;
	const float mu = hs->u1 + dy21 * (hs->u3 - hs->u1) * inv_dy31;
	const float mv = hs->v1 + dy21 * (hs->v3 - hs->v1) * inv_dy31;
	const float mw = hs->w1 + dy21 * (hs->w3 - hs->w1) * inv_dy31;

	if (hs->sc < hs->sb)
	{
		hs->dudx = slope(mu, mx, hs->u2, hs->p2->x, UV_SHIFT);
		hs->dudy = slope(hs->u1, hs->p1->y, hs->u3, hs->p3->y, UV_SHIFT);
		hs->dvdx = slope(mv, mx, hs->v2, hs->p2->x, UV_SHIFT);
		hs->dvdy = slope(hs->v1, hs->p1->y, hs->v3, hs->p3->y, UV_SHIFT);
		hs->dwdx = slope(mw, mx, hs->w2, hs->p2->x, W_SHIFT);
		hs->dwdy = slope(hs->w1, hs->p1->y, hs->w3, hs->p3->y, W_SHIFT);
	}
	else
	{
		hs->dudx = slope(hs->u2, hs->p2->x, mu, mx, UV_SHIFT);
		hs->dudy = slope(hs->u1, hs->p1->y, hs->u2, hs->p2->y, UV_SHIFT);
		hs->dvdx = slope(hs->v2, hs->p2->x, mv, mx, UV_SHIFT);
		hs->dvdy = slope(hs->v1, hs->p1->y, hs->v2, hs->p2->y, UV_SHIFT);
		hs->dwdx = slope(hs->w2, hs->p2->x, mw, mx, W_SHIFT);
		hs->dwdy = slope(hs->w1, hs->p1->y, hs->w2, hs->p2->y, W_SHIFT);
	}

	hs->urow = (int)(hs->u1 * (1 << UV_SHIFT));
	hs->vrow = (int)(hs->v1 * (1 << UV_SHIFT));
	hs->wrow = (int)(hs->w1 * (1 << W_SHIFT));
	return true;
}

static void switch_to_part_23(raster_scanline_t *hs)
{
	assert(hs->part_12);
	hs->part_12 = false;

	int new_dx = slope(hs->p2->x, hs->p2->y, hs->p3->x, hs->p3->y, 16);
	if (hs->sb < hs->sc)
	{
		hs->dx1 = new_dx;
		hs->dudy = slope(hs->u2, hs->p2->y, hs->u3, hs->p3->y, UV_SHIFT);
		hs->dvdy = slope(hs->v2, hs->p2->y, hs->v3, hs->p3->y, UV_SHIFT);
		hs->x1row = (int)(hs->p2->x * (1 << 16));
		hs->urow = (int)(hs->u2 * (1 << UV_SHIFT));
		hs->vrow = (int)(hs->v2 * (1 << UV_SHIFT));
		hs->dwdy = slope(hs->w2, hs->p2->y, hs->w3, hs->p3->y, W_SHIFT);
		hs->wrow = (int)(hs->w2 * (1 << W_SHIFT));
	}
	else
	{
		hs->x2row = (int)(hs->p2->x * (1 << 16));
		hs->dx2 = new_dx;
	}
}

FORCE_INLINE bool raster_scanline_y_continue(raster_scanline_t* hs)
{
	if (hs->part_12 && hs->y == hs->y2)
	{
		switch_to_part_23(hs);
	}
	const int cur_y_end = hs->part_12 ? hs->y2 : hs->y3;
	if (cur_y_end < 0)
	{
		assert(hs->part_12);
		int dy = cur_y_end - hs->y;
		hs->x1row += dy * hs->dx1;
		hs->x2row += dy * hs->dx2;
		hs->urow += dy * hs->dudy;
		hs->vrow += dy * hs->dvdy;
		hs->wrow += dy * hs->dwdy;
		hs->y = cur_y_end;
		switch_to_part_23(hs);
	}
	if (hs->y < 0)
	{
		hs->x1row -= hs->y * hs->dx1;
		hs->x2row -= hs->y * hs->dx2;
		hs->urow -= hs->y * hs->dudy;
		hs->vrow -= hs->y * hs->dvdy;
		hs->wrow -= hs->y * hs->dwdy;
		hs->y = 0;
	}

	return hs->y < hs->y3;
}

FORCE_INLINE void raster_scanline_y_step(raster_scanline_t* hs)
{
	hs->x1row += hs->dx1;
	hs->x2row += hs->dx2;
	hs->urow += hs->dudy;
	hs->vrow += hs->dvdy;
	hs->wrow += hs->dwdy;
	hs->y++;
}

FORCE_INLINE void raster_scanline_x_begin(raster_scanline_t* hs)
{
	const int x1 = hs->x1row >> 16;
	hs->endx = (hs->x2row >> 16) + 1;
	hs->x = x1;
	hs->u = hs->urow;
	hs->v = hs->vrow;
	hs->w = hs->wrow;

	hs->endx = min2(hs->endx, SCREEN_X);

	if (hs->x < 0)
	{
		hs->u -= hs->x * hs->dudx;
		hs->v -= hs->x * hs->dvdx;
		hs->w -= hs->x * hs->dwdx;
		hs->x = 0;
	}
}

FORCE_INLINE bool raster_scanline_x_continue(raster_scanline_t* hs)
{
	return hs->x < hs->endx;
}

FORCE_INLINE void raster_scanline_x_step(raster_scanline_t* hs)
{
	hs->u += hs->dudx;
	hs->v += hs->dvdx;
	hs->w += hs->dwdx;
	hs->x++;
}

// --------------------------------------------------------------------------

static void draw_triangle_pattern_scanline(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6], const uint8_t *pattern)
{
	raster_scanline_t hs;
	if (!raster_scanline_begin(&hs, p1, p2, p3, uvs, 64.0f))
		return;

	for (; raster_scanline_y_continue(&hs); raster_scanline_y_step(&hs))
	{
		uint8_t* row = bitmap + hs.y * rowstride;
		raster_scanline_x_begin(&hs);
		if (!raster_scanline_x_continue(&hs))
			continue;

		uint8_t pat = pattern[hs.y % 8];
		uint32_t color = (pat << 24) | (pat << 16) | (pat << 8) | pat;

		// write out pixels 32 at a time
		int startbit = hs.x % 32;
		uint32_t startmask = swap((1 << (32 - startbit)) - 1);
		int endbit = hs.endx % 32;
		uint32_t endmask = swap(((1 << endbit) - 1) << (32 - endbit));

		int col = hs.x / 32;
		uint32_t* p = (uint32_t*)row + col;

		if (col == hs.endx / 32)
		{
			uint32_t mask = 0;
			if (startbit > 0 && endbit > 0)
				mask = startmask & endmask;
			else if (startbit > 0)
				mask = startmask;
			else if (endbit > 0)
				mask = endmask;
			_drawMaskPattern(p, mask, color);
		}
		else
		{
			int x = hs.x;
			if (startbit > 0)
			{
				_drawMaskPattern(p++, startmask, color);
				x += (32 - startbit);
			}
			while (x + 32 <= hs.endx)
			{
				_drawMaskPattern(p++, 0xffffffff, color);
				x += 32;
			}
			if (endbit > 0)
			{
				_drawMaskPattern(p, endmask, color);
			}
		}
	}
}

static void draw_triangle_bluenoise_scanline(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6], const uint8_t tri_color)
{
	raster_scanline_t hs;
	if (!raster_scanline_begin(&hs, p1, p2, p3, uvs, 64.0f))
		return;

	const uint8_t* blue_noise_buffer = get_blue_noise_buffer();

	for (; raster_scanline_y_continue(&hs); raster_scanline_y_step(&hs))
	{
		raster_scanline_x_begin(&hs);
		if (!raster_scanline_x_continue(&hs))
			continue;

		const uint8_t* noise_row = blue_noise_buffer + hs.y * SCREEN_X;
		uint8_t* row = bitmap + hs.y * rowstride;

		int x = hs.x, endx = hs.endx;

		// write out pixels 32 at a time
		uint32_t mask = 0;
		uint32_t* p = (uint32_t*)row + hs.x / 32;
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
}

static void draw_triangle_checker_scanline(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6])
{
	raster_scanline_t hs;
	if (!raster_scanline_begin(&hs, p1, p2, p3, uvs, 64.0f))
		return;

	for (; raster_scanline_y_continue(&hs); raster_scanline_y_step(&hs))
	{
		uint8_t* row = bitmap + hs.y * rowstride;
		raster_scanline_x_begin(&hs);

		// write out pixels 32 at a time
		uint32_t mask = 0;
		uint32_t* p = (uint32_t*)row + hs.x / 32; // 1.0%
		uint32_t color = 0;

		while(raster_scanline_x_continue(&hs)) // 7.1%
		{
			mask |= 0x80000000u >> (hs.x & 31); // 8%

			int uu = hs.u * 5; // 0.4%
			int vv = hs.v * 5;
			// |1 to prevent floating point division error
			int divisor = (hs.w >> max2(0, W_SHIFT - UV_SHIFT)) | 1;
			uint16_t ui = (uu / divisor) >> max2(0, UV_SHIFT - W_SHIFT);
			uint16_t vi = (vv / divisor) >> max2(0, UV_SHIFT - W_SHIFT);
			ui &= 63;
			vi &= 63;
			bool checker = (ui > 31) != (vi > 31);

			uint32_t texpix = checker ? 0 : 0x80;
			color |= (texpix << 24) >> (hs.x % 32);

			//int bit_mask = 1 << (7 - (hs.x & 7));
			//if (checker)
			//	row[hs.x / 8] &= ~bit_mask;
			//else
			//	row[hs.x / 8] |= bit_mask;

			raster_scanline_x_step(&hs);

			if (hs.x % 32 == 0) // 5.6%
			{
				_drawMaskPattern(p++, swap(mask), swap(color));
				mask = 0;
				color = 0;
			}
		}

		_drawMaskPattern(p, swap(mask), swap(color));
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

	float dxu, dxv, dyu, dyv;

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
	state->dyu = state->u_01 - state->u_00;
	state->dyv = state->v_01 - state->v_00;
	state->dxu = state->u_10 - state->u_00;
	state->dxv = state->v_10 - state->v_00;
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
#define DITHER_RES_MASK (DITHER_RES-1)
#define DITHER_DOTS_PER_SIDE (4) // Note: upstream used (DITHER_RES/16)
#define DITHER_SLICES (DITHER_DOTS_PER_SIDE * DITHER_DOTS_PER_SIDE)

#define DITHER_SCALE (6.0f)
#define DITHER_SCALE_EXP (64.0f) // exp2f(DITHER_SCALE)

FORCE_INLINE bool do_pixel_sample(const float u, const float v, const int patternScaleLevel_i, const int subLayer_offset, const int compare, const int x, const int y)
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
	// Get the pattern value relative to the threshold, scale it
	// according to the contrast, and add the base value.
	// Note: the terms are slightly reassociated compared to upstream, since
	// we only need to threshold the result.

#if 0
	uint8_t* debug_output = plat_gfx_get_debug_frame();
	if (debug_output) {
		float3 dbg = (float3){ 0,0,0 };

		// UV
		dbg.x = fract(u), dbg.y = fract(v);

		// fractal UV
		//dbg.x = fract(uu), dbg.y = fract(vv);

		// pattern
		dbg.x = dbg.y = dbg.z = pattern / 255.0f;

		// bw
		//dbg.x = dbg.y = dbg.z = check < 0.0f ? 0.0f : 1.0f;

		const int offset = (y * SCREEN_X + x) * 4 + 0;
		debug_output[offset + 0] = (uint8_t)(saturate(dbg.x) * 255.0f);
		debug_output[offset + 1] = (uint8_t)(saturate(dbg.y) * 255.0f);
		debug_output[offset + 2] = (uint8_t)(saturate(dbg.z) * 255.0f);
		debug_output[offset + 3] = 255;
	}
#endif
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
			float spacing = (fabsf(state.dxu) + fabsf(state.dxv) + fabsf(state.dyu) + fabsf(state.dyv)) * 0.25f;
			spacing *= spacingMul;

			int patternScaleLevel_i;
			const int subLayer_offset = get_dither3d_level_fraction(spacing, &patternScaleLevel_i);

			// Note: accumulating pixel outputs/masks into 32 bit words and writing them out
			// once they are done (or trailing after the loop) seems to be slightly slower on Playdate,
			// than simply operating on bytes in a naive way.
			const int byte_idx = state.x / 8;
			const int mask_0 = 1 << (7 - ((state.x + 0) & 7));
			const int mask_1 = 1 << (7 - ((state.x + 1) & 7));

			if (state.inside_00)
			{
				bool check = do_pixel_sample(state.u_00, state.v_00, patternScaleLevel_i, subLayer_offset, compare_val, state.x, state.y);
				if (check)
					output_0[byte_idx] &= ~mask_0;
				else
					output_0[byte_idx] |= mask_0;
			}
			if (state.inside_01)
			{
				bool check = do_pixel_sample(state.u_01, state.v_01, patternScaleLevel_i, subLayer_offset, compare_val, state.x + 1, state.y);
				if (check)
					output_0[byte_idx] &= ~mask_1;
				else
					output_0[byte_idx] |= mask_1;
			}
			if (state.inside_10)
			{
				bool check = do_pixel_sample(state.u_10, state.v_10, patternScaleLevel_i, subLayer_offset, compare_val, state.x, state.y + 1);
				if (check)
					output_1[byte_idx] &= ~mask_0;
				else
					output_1[byte_idx] |= mask_0;
			}
			if (state.inside_11)
			{
				bool check = do_pixel_sample(state.u_11, state.v_11, patternScaleLevel_i, subLayer_offset, compare_val, state.x + 1, state.y + 1);
				if (check)
					output_1[byte_idx] &= ~mask_1;
				else
					output_1[byte_idx] |= mask_1;
			}
		}
	}
}

// --------------------------------------------------------------------------
// Simple pure black/white checkerboard with halfspace rasterizer

FORCE_INLINE bool do_checker(float u, float v)
{
	u = fract(u * 5.0f);
	v = fract(v * 5.0f);
	return (u > 0.5f) != (v > 0.5f);
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
				bool check = do_checker(state.u_00, state.v_00);
				if (check)
					output_0[byte_idx] &= ~mask_0;
				else
					output_0[byte_idx] |= mask_0;
			}
			if (state.inside_01)
			{
				bool check = do_checker(state.u_01, state.v_01);
				if (check)
					output_0[byte_idx] &= ~mask_1;
				else
					output_0[byte_idx] |= mask_1;
			}
			if (state.inside_10)
			{
				bool check = do_checker(state.u_10, state.v_10);
				if (check)
					output_1[byte_idx] &= ~mask_0;
				else
					output_1[byte_idx] |= mask_0;
			}
			if (state.inside_11)
			{
				bool check = do_checker(state.u_11, state.v_11);
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

static void drawShapeFace(const Scene* scene, uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float3* normal, const Mesh* mesh, int tri_index, enum DrawStyle style, bool wire)
{
	// If any vertex is behind the camera, skip it
	if (p1->z <= 0 || p2->z <= 0 || p3->z <= 0)
		return;

	float x1 = p1->x;
	float y1 = p1->y;
	float x2 = p2->x;
	float y2 = p2->y;
	float x3 = p3->x;
	float y3 = p3->y;

	// quick bounds check
	if ((x1 < 0 && x2 < 0 && x3 < 0) ||
		(x1 >= SCREEN_X && x2 >= SCREEN_X && x3 >= SCREEN_X) ||
		(y1 < 0 && y2 < 0 && y3 < 0) ||
		(y1 >= SCREEN_Y && y2 >= SCREEN_Y && y3 >= SCREEN_Y))
		return;

	// only render front side of faces via winding order
	float dx21 = x2 - x1;
	float dy31 = y3 - y1;
	float dx31 = x3 - x1;
	float dy21 = y2 - y1;
	float d = dx21 * dy31 - dy21 * dx31;
	if (d >= 0)
		return;

	//float kSmallPx = 8.0f;
	//if (fabsf(dx21) < kSmallPx && fabsf(dy31) < kSmallPx && fabsf(dx31) < kSmallPx && fabsf(dy21) < kSmallPx)
	//	return;

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
		draw_triangle_pattern_scanline(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6, pattern);
	}
	else if (style == Draw_Bluenoise)
	{
		// draw into byte buffer for blue noise thresholding
		int col = (int)(v * 255.0f);
		if (col < 0) col = 0;
		if (col > 255) col = 255;

		draw_triangle_bluenoise_scanline(bitmap, rowstride, p1, p2, p3, mesh->uvs + tri_index * 6, col);
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

	// sort faces by z
	s_facesort_instance = scene;
	qsort(scene->tmp_order_table, mesh->tri_count, sizeof(scene->tmp_order_table[0]), compareFaceZ);

	// draw faces
	for (int i = 0; i < mesh->tri_count; ++i)
	{
		uint16_t fi = scene->tmp_order_table[i];
		uint16_t idx0 = mesh->tris[fi * 3 + 0];
		uint16_t idx1 = mesh->tris[fi * 3 + 1];
		uint16_t idx2 = mesh->tris[fi * 3 + 2];
		drawShapeFace(scene, buffer, rowstride, &scene->tmp_points[idx0], &scene->tmp_points[idx1], &scene->tmp_points[idx2], &scene->tmp_face_normals[fi], mesh, fi, style, wire);
	}
}
