#pragma once

#include "mathlib.h"

typedef struct Mesh {
	int vertex_count;
	const float3* vertices;
	int tri_count;
	const uint16_t* tris;
	const float* uvs; // 6x float per triangle
} Mesh;

typedef struct Scene
{
	xform camera;
	float3 light;

	// location of the Z vanishing point on the screen. (0,0) is top left corner, (1,1) is bottom right. Defaults to (0.5,0.5)
	float centerx;
	float centery;

	// display scaling factor. Default is 120, so that the screen extents are (-1.66,1.66)x(-1,1)
	float scale;

	int tmp_points_cap;
	int tmp_faces_cap;
	float3* tmp_points;
	float3* tmp_face_normals;
	uint16_t* tmp_order_table;
	float* tmp_dist_table;
} Scene;

typedef struct MeshInstance {
	const Mesh* mesh;
	xform tr;
} MeshInstance;

void drawLine(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, int thick, const uint8_t pattern[8]);
void drawShapeFace(const Scene* scene, uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float3* normal, const Mesh* mesh, int tri_index, enum DrawStyle style, bool wire);

void scene_init(Scene* scene);
void scene_shutdown(Scene* scene);
void scene_setCamera(Scene* scene, float3 origin, float3 lookAt, float scale, float3 up);
void scene_setLight(Scene* scene, float3 light);
void scene_setCenter(Scene* scene, float x, float y);

enum DrawStyle {
	// 15ms ditherpattern from mini3d
	// 19ms generic scanline rasterizer, that fills scanlines with dither pattern (but still has per-row overhead of UVs etc.)
	Draw_Pattern,
	// 31ms bluenoise dither, rasterizer loop from mini3
	// 33ms bluenoise dither, generic scanline rasterizer loop (more per-row overhead)
	// 24ms no bluenoise screen buffer; use the bluenoise thresholding directly in rasterizer
	Draw_Bluenoise,

	// 35ms device, 0.46mc PC
	Draw_Checker_Scanline,

	// 42ms device, fl_fl_div, write 32 pixels
	Draw_Checker_Hecker,
	Draw_Checker_HeckerSub,

	// 53ms device halfspace
	Draw_Checker_Halfspace,

	Draw_Dither3D_Scanline, // wrong/incorrect right now
	// 102ms device, fl_fl_div, perspective correct derivatives, write each bit
	Draw_Dither3D_Hecker,

	// 605ms device, 5.8ms PC
	// 316ms device, 4.2ms PC adjust_float_exp
	// 246ms device, 2.7ms PC bit manip for patternScaleLevel and f
	// 149ms device, 1.4ms 2x2 raster
	// 137ms device, 1.3ms simplification to contrast maths
	// 123ms device, 1.4ms move final sampling inside visibility check for each 2x2 pixel
	// 107ms device, 1.3ms simplify spacing & contrast calc (behavior change!)
	// 83ms device, reduce dither texture resolution 2x on XY
	// 80ms device, enable -ffast-math
	// 76ms device, perspective correct UVs at 2x horizontal step
	// 74ms device, skip rest of row once we get out of triangle
	// 65ms device, simplify sampled pattern comparison
	// 63ms device, small micro opts
	Draw_Dither3D_Halfspace,
	Draw_Count
};

void scene_drawMesh(Scene* scene, uint8_t* buffer, int rowstride, const Mesh* mesh, const xform* matrix, enum DrawStyle style, bool wire);
