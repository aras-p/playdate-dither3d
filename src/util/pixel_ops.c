// SPDX-License-Identifier: Unlicense

#include "pixel_ops.h"
#include "image_loader.h"

#include "../globals.h"
#include "../mathlib.h"

#include "../platform.h"

#include <string.h>
#include <stdlib.h>

static uint8_t* s_blue_noise;
uint8_t* s_dither4x4_r;
uint8_t* s_dither4x4_g;

uint8_t g_screen_buffer[SCREEN_X * SCREEN_Y];
uint8_t g_screen_buffer_2x2sml[SCREEN_X/2 * SCREEN_Y/2];

void clear_screen_buffers()
{
	memset(g_screen_buffer, 0xFF, sizeof(g_screen_buffer));
	memset(g_screen_buffer_2x2sml, 0xFF, sizeof(g_screen_buffer_2x2sml));
}


// 2x2 pixel block ordered dither matrix.
// 0 3
// 2 1
int g_order_pattern_2x2[4][2] = {
	{1, 0},
	{0, 2},
	{2, 0},
	{0, 1},
};

// 3x2 pixel block ordered dither matrix.
// 0 2 4
// 3 5 1
int g_order_pattern_3x2[6][2] = {
	{1, 0},
	{0, 3},
	{2, 0},
	{0, 1},
	{3, 0},
	{0, 2},
};

// 4x2 pixel block ordered dither matrix.
// 0 4 2 6
// 3 7 1 5
int g_order_pattern_4x2[8][2] = {
	{1, 0},
	{0, 3},
	{3, 0},
	{0, 1},
	{2, 0},
	{0, 4},
	{4, 0},
	{0, 2},
};

// 4x3 pixel block ordered dither matrix.
// 0 9 6 3
// 7 4 8 B
// 2 A 1 5
int g_order_pattern_4x3[12][3] = {
	{1, 0, 0},
	{0, 0, 3},
	{0, 0, 1},
	{4, 0, 0},
	{0, 2, 0},
	{0, 0, 4},
	{3, 0, 0},
	{0, 1, 0},
	{0, 3, 0},
	{2, 0, 0},
	{0, 0, 2},
	{0, 4, 0},
};

// 4x4 pixel block ordered dither matrix.
//  0 12  3 15
//  8  4 11  7
//  2 14  1 13
// 10  6  9  5
int g_order_pattern_4x4[16][4] = {
	{1, 0, 0, 0},
	{0, 0, 3, 0},
	{0, 0, 1, 0},
	{3, 0, 0, 0},
	{0, 2, 0, 0},
	{0, 0, 0, 4},
	{0, 0, 0, 2},
	{0, 4, 0, 0},
	{0, 1, 0, 0},
	{0, 0, 0, 3},
	{0, 0, 0, 1},
	{0, 3, 0, 0},
	{2, 0, 0, 0},
	{0, 0, 4, 0},
	{0, 0, 2, 0},
	{4, 0, 0, 0},
};

void init_pixel_ops()
{
	{
		int ww, hh;
		s_blue_noise = read_tga_file_grayscale("BlueNoise.tga", &ww, &hh);
		if (ww != SCREEN_X || hh != SCREEN_Y) {
			plat_free(s_blue_noise);
			s_blue_noise = NULL;
		}
	}
	{
		int ww, hh;
		s_dither4x4_r = read_tga_file_grayscale("Dither3D_4x4-sm-r.tga", &ww, &hh);
		if (ww != 64/2 || hh != 1024/2) {
			plat_free(s_dither4x4_r);
			s_dither4x4_r = NULL;
		}
	}
	{
		int ww, hh;
		s_dither4x4_g = read_tga_file_grayscale("Dither3D_4x4-g.tga", &ww, &hh);
		if (ww != 64 || hh != 1) {
			plat_free(s_dither4x4_g);
			s_dither4x4_g = NULL;
		}
	}

	memset(g_screen_buffer, 0xFF, sizeof(g_screen_buffer));
	memset(g_screen_buffer_2x2sml, 0xFF, sizeof(g_screen_buffer_2x2sml));
}


void draw_dithered_scanline(const uint8_t* values, int y, int bias, uint8_t* framebuffer)
{
	const uint32_t* values32 = (const uint32_t*)values;
	uint8_t* row = framebuffer + y * SCREEN_STRIDE_BYTES;
	const uint32_t* noise_row = (const uint32_t*)(s_blue_noise + y * SCREEN_X);
	for (int bx = 0; bx < SCREEN_X / 8; ++bx) {
		uint32_t noise03 = *noise_row++;
		uint32_t noise47 = *noise_row++;
		uint32_t values03 = *values32++;
		uint32_t values47 = *values32++;
		uint8_t pixbyte = 0;
		if ((values03 & 0xFF) > (noise03 & 0xFF) + bias) pixbyte |= 128;
		if ((values03 & 0xFF00) > (noise03 & 0xFF00) + bias) pixbyte |= 64;
		if ((values03 & 0xFF0000) > (noise03 & 0xFF0000) + bias) pixbyte |= 32;
		if ((values03 & 0xFF000000) > (noise03 & 0xFF000000) + bias) pixbyte |= 16;
		if ((values47 & 0xFF) > (noise47 & 0xFF) + bias) pixbyte |= 8;
		if ((values47 & 0xFF00) > (noise47 & 0xFF00) + bias) pixbyte |= 4;
		if ((values47 & 0xFF0000) > (noise47 & 0xFF0000) + bias) pixbyte |= 2;
		if ((values47 & 0xFF000000) > (noise47 & 0xFF000000) + bias) pixbyte |= 1;
		*row++ = pixbyte;
	}
}

void draw_dithered_screen(uint8_t* framebuffer, int bias)
{
	const uint8_t* src = g_screen_buffer;
	for (int y = 0; y < SCREEN_Y; ++y)
	{
		draw_dithered_scanline(src, y, bias, framebuffer);
		src += SCREEN_X;
	}
}

void draw_dithered_screen_2x2(uint8_t* framebuffer, int filter)
{
	uint8_t rowvalues[SCREEN_X];
	if (filter == 0)
	{
		// use same source value for each 2x2 block
		int src_idx = 0;
		for (int y = 0; y < SCREEN_Y / 2; ++y) {
			for (int x = 0; x < SCREEN_X / 2; ++x, ++src_idx) {
				uint8_t val = g_screen_buffer[src_idx];
				rowvalues[x * 2 + 0] = val;
				rowvalues[x * 2 + 1] = val;
			}
			draw_dithered_scanline(rowvalues, y * 2 + 0, 0, framebuffer);
			draw_dithered_scanline(rowvalues, y * 2 + 1, 0, framebuffer);
		}
	}
	else if (filter == 1)
	{
		// filter values horizontally: bottom to top, right to left.
		// this will expand image to be full width, but half height.
		for (int y = SCREEN_Y/2 - 1; y >= 0; --y) {
			int src_idx = y * SCREEN_X / 2;
			int dst_idx = y * SCREEN_X;
			for (int x = SCREEN_X / 2 - 1; x >= 0; --x) {
				int val_prev = g_screen_buffer[src_idx + (x > 0 ? x - 1 : 0)];
				int val_curr = g_screen_buffer[src_idx + x];
				g_screen_buffer[dst_idx + x * 2 + 0] = (val_prev + val_curr) >> 1;
				g_screen_buffer[dst_idx + x * 2 + 1] = val_curr;
			}
		}
		// draw scanlines, filtering odd ones
		for (int y = 0; y < SCREEN_Y / 2; ++y)
		{
			const uint8_t* row = &g_screen_buffer[y * SCREEN_X];
			// draw raw scanline
			draw_dithered_scanline(row, y * 2 + 0, 0, framebuffer);
			int next_y = y * 2 + 2;
			if (next_y >= SCREEN_Y)
			{
				// nothing to filter with below, just draw previous
				draw_dithered_scanline(row, y * 2 + 1, 0, framebuffer);
			}
			else
			{
				// compute filtered scanline
				for (int x = 0; x < SCREEN_X; ++x)
				{
					rowvalues[x] = ((int)row[x] + (int)row[x + SCREEN_X]) >> 1;
				}
				draw_dithered_scanline(rowvalues, y * 2 + 1, 0, framebuffer);
			}
		}
	}
}

// DDA line drawing algorithm, using 16.16 fixed point
void draw_line(uint8_t* framebuffer, int width, int height, int x1, int y1, int x2, int y2, uint8_t color)
{
	int dx = x2 - x1;
	int dy = y2 - y1;
	int abs_dx = abs(dx);
	int abs_dy = abs(dy);
	int steps = max2(abs_dx, abs_dy);
	if (steps == 0)
		return;

	int xstep_fx = (dx << 16) / steps;
	int ystep_fx = (dy << 16) / steps;

	int x_fx = x1 << 16;
	int y_fx = y1 << 16;

	for (int i = 0; i <= steps; ++i)
	{
		int ix = x_fx >> 16;
		int iy = y_fx >> 16;
		if (ix >= 0 && ix < width && iy >= 0 && iy < height)
			framebuffer[iy * width + ix] = color;
		x_fx += xstep_fx;
		y_fx += ystep_fx;
	}
}
