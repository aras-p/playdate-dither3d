#include "draw_style.h"
#include "rasterizer_halfspace.h"
#include "rasterizer_scanline.h"
#include "util/pixel_ops.h"

FORCE_INLINE void _drawMaskPattern(uint32_t* p, uint32_t mask, uint32_t color)
{
	if ( mask == 0xffffffff )
		*p = color;
	else
		*p = (*p & ~mask) | (color & mask);
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

void draw_tri_pattern_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t* pattern)
{
	raster_scanline_simple_t tri;
	raster_scanline_simple_begin(&tri, p0, p1, p2);

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

void draw_tri_bluenoise_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t tri_color)
{
	raster_scanline_simple_t tri;
	raster_scanline_simple_begin(&tri, p0, p1, p2);

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


void draw_tri_dither3d_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6], const uint8_t brightness)
{
	raster_halfspace_t hs;
	raster_halfspace_begin(&hs, p1, p2, p3, uvs);

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
			uint32_t texpix = check ? 0x0 : 0x80000000u;
			color |= texpix >> (X & 31);


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
			uint32_t texpix = check ? 0x0 : 0x80000000u;
			color |= texpix >> (X & 31);

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

void draw_tri_dither3d_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const uint8_t brightness)
{
	raster_scanline_uv_t tri;
	raster_scanline_uv_begin(&tri, p0, p1, p2, uvs, 1.0f);

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

void draw_tri_checker_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6])
{
	raster_scanline_uv_t tri;
	raster_scanline_uv_begin(&tri, p0, p1, p2, uvs, 5.0f);

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

void draw_tri_checker_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6])
{
	raster_halfspace_t hs;
	raster_halfspace_begin(&hs, p1, p2, p3, uvs);

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
