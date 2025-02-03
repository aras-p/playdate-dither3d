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
void fillTriangle(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const uint8_t pattern[8]);

void scene_init(Scene* scene);
void scene_shutdown(Scene* scene);
void scene_setCamera(Scene* scene, float3 origin, float3 lookAt, float scale, float3 up);
void scene_setLight(Scene* scene, float3 light);
void scene_setCenter(Scene* scene, float x, float y);

enum DrawStyle {
	// 45fps
	Draw_Pattern = 0,
	// 26fps (32 no wire)
	Draw_Bluenoise,

	// 26fps simple (as bluenoise, just with halfspace raster)
	// 17fps (20fps no wire) checker, no perspective correction
	// 15fps (17fps no wire) checker, with perspective correction
	// 18fps (20fps no wire) checker, with perspective, no blue noise just b/w
	Draw_Checker_Halfspace,
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
	Draw_Dither3D_Halfspace,
	Draw_Count
};

void scene_drawMesh(Scene* scene, uint8_t* buffer, int rowstride, const Mesh* mesh, const xform* matrix, enum DrawStyle style, bool wire);
