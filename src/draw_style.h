#pragma once

#include "mathlib.h"

void draw_tri_pattern_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t* pattern);
void draw_tri_bluenoise_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const uint8_t tri_color);

void draw_tri_checker_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6]);
void draw_tri_checker_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6]);

void draw_tri_dither3d_scanline(uint8_t* bitmap, int rowstride, const float3* p0, const float3* p1, const float3* p2, const float uvs[6], const uint8_t brightness);
void draw_tri_dither3d_halfspace(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float uvs[6], const uint8_t brightness);
