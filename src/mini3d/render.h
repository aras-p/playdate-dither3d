#pragma once

#include "../mathlib.h"

typedef struct Mesh {
	int vertex_count;
	const float3* vertices;
	int tri_count;
	const uint16_t* tris;
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


void drawLine(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, int thick, const uint8_t pattern[8]);
void fillTriangle(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const uint8_t pattern[8]);

void scene_init(Scene* scene);
void scene_shutdown(Scene* scene);
void scene_setCamera(Scene* scene, float3 origin, float3 lookAt, float scale, float3 up);
void scene_setLight(Scene* scene, float3 light);
void scene_setCenter(Scene* scene, float x, float y);

void scene_drawMesh(Scene* scene, uint8_t* buffer, int rowstride, const Mesh* mesh, const xform* matrix);
