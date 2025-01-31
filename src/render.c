// Parts based on mini-3d example from Playdate SDK:
// Created by Dave Hayden on 10/20/15.
// Copyright © 2015 Panic, Inc. All rights reserved.

#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "render.h"
#include "platform.h"
#include "util/pixel_ops.h"

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

static inline int32_t slope(float x1, float y1, float x2, float y2)
{
	float dx = x2-x1;
	float dy = y2-y1;
	
	if ( dy < 1 )
		return (int32_t)(dx * (1<<16));
	else
		return (int32_t)(dx / dy * (1<<16));
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
	int32_t dx = slope(p1->x, p1->y, p2->x, p2->y);
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

static void fillRange(uint8_t* bitmap, int rowstride, int y, int endy, int32_t* x1p, int32_t dx1, int32_t* x2p, int32_t dx2, const uint8_t pattern[8])
{
	int32_t x1 = *x1p, x2 = *x2p;
	
	if ( endy < 0 )
	{
		int dy = endy - y;
		*x1p = x1 + dy * dx1;
		*x2p = x2 + dy * dx2;
		return;
	}
	
	if ( y < 0 )
	{
		x1 += -y * dx1;
		x2 += -y * dx2;
		y = 0;
	}
	
	while ( y < endy )
	{
		uint8_t p = pattern[y%8];
		uint32_t color = (p<<24) | (p<<16) | (p<<8) | p;
		
		drawFragment((uint32_t*)&bitmap[y*rowstride], (x1>>16), (x2>>16)+1, color);
		
		x1 += dx1;
		x2 += dx2;
		++y;
	}
	
	*x1p = x1;
	*x2p = x2;
}

static inline void sortTri(const float3** p1, const float3** p2, const float3** p3)
{
	float y1 = (*p1)->y, y2 = (*p2)->y, y3 = (*p3)->y;
	
	if ( y1 <= y2 && y1 < y3 )
	{
		if ( y3 < y2 ) // 1,3,2
		{
			const float3* tmp = *p2;
			*p2 = *p3;
			*p3 = tmp;
		}
	}
	else if ( y2 < y1 && y2 < y3 )
	{
		const float3* tmp = *p1;
		*p1 = *p2;

		if ( y3 < y1 ) // 2,3,1
		{
			*p2 = *p3;
			*p3 = tmp;
		}
		else // 2,1,3
			*p2 = tmp;
	}
	else
	{
		const float3* tmp = *p1;
		*p1 = *p3;
		
		if ( y1 < y2 ) // 3,1,2
		{
			*p3 = *p2;
			*p2 = tmp;
		}
		else // 3,2,1
			*p3 = tmp;
	}

}

void fillTriangle(uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const uint8_t pattern[8])
{
	// sort by y coord
	
	sortTri(&p1, &p2, &p3);
	
	int endy = min2(SCREEN_Y, (int)p3->y);
	
	if ( p1->y > SCREEN_Y || endy < 0 )
		return;

	int32_t x1 = (int32_t)(p1->x * (1<<16));
	int32_t x2 = x1;
	
	int32_t sb = slope(p1->x, p1->y, p2->x, p2->y);
	int32_t sc = slope(p1->x, p1->y, p3->x, p3->y);

	int32_t dx1 = min2(sb, sc);
	int32_t dx2 = max2(sb, sc);
	
	fillRange(bitmap, rowstride, (int)p1->y, min2(SCREEN_Y, (int)p2->y), &x1, dx1, &x2, dx2, pattern);
	
	int dx = slope(p2->x, p2->y, p3->x, p3->y);
	
	if ( sb < sc )
	{
		x1 = (int32_t)(p2->x * (1<<16));
		fillRange(bitmap, rowstride, (int)p2->y, endy, &x1, dx, &x2, dx2, pattern);
	}
	else
	{
		x2 = (int32_t)(p2->x * (1<<16));
		fillRange(bitmap, rowstride, (int)p2->y, endy, &x1, dx1, &x2, dx, pattern);
	}
}

static void draw_fragment(uint8_t* row, int x1, int x2, const uint8_t color)
{
	if (x2 < 0 || x1 >= SCREEN_X)
		return;

	if (x1 < 0)
		x1 = 0;

	if (x2 > SCREEN_X)
		x2 = SCREEN_X;

	if (x1 > x2)
		return;

	memset(row + x1, color, x2 - x1);
}

static void fill_range(int y, int endy, int32_t* x1p, int32_t dx1, int32_t* x2p, int32_t dx2, const uint8_t color)
{
	int32_t x1 = *x1p, x2 = *x2p;

	if (endy < 0)
	{
		int dy = endy - y;
		*x1p = x1 + dy * dx1;
		*x2p = x2 + dy * dx2;
		return;
	}

	if (y < 0)
	{
		x1 += -y * dx1;
		x2 += -y * dx2;
		y = 0;
	}

	while (y < endy)
	{
		draw_fragment(g_screen_buffer + y * SCREEN_X, (x1 >> 16), (x2 >> 16) + 1, color);

		x1 += dx1;
		x2 += dx2;
		++y;
	}

	*x1p = x1;
	*x2p = x2;
}

void draw_triangle(const float3* p1, const float3* p2, const float3* p3, const uint8_t col)
{
	// sort by y coord
	sortTri(&p1, &p2, &p3);

	int endy = min2(SCREEN_Y, (int)p3->y);

	if (p1->y > SCREEN_Y || endy < 0)
		return;

	int32_t x1 = (int32_t)(p1->x * (1 << 16));
	int32_t x2 = x1;

	int32_t sb = slope(p1->x, p1->y, p2->x, p2->y);
	int32_t sc = slope(p1->x, p1->y, p3->x, p3->y);

	int32_t dx1 = min2(sb, sc);
	int32_t dx2 = max2(sb, sc);

	fill_range((int)p1->y, min2(SCREEN_Y, (int)p2->y), &x1, dx1, &x2, dx2, col);

	int dx = slope(p2->x, p2->y, p3->x, p3->y);

	if (sb < sc)
	{
		x1 = (int32_t)(p2->x * (1 << 16));
		fill_range((int)p2->y, endy, &x1, dx, &x2, dx2, col);
	}
	else
	{
		x2 = (int32_t)(p2->x * (1 << 16));
		fill_range((int)p2->y, endy, &x1, dx1, &x2, dx, col);
	}
}

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
	scene->light = light;
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

static void drawShapeFace(const Scene* scene, uint8_t* bitmap, int rowstride, const float3* p1, const float3* p2, const float3* p3, const float3* normal, enum DrawStyle style, bool wire)
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
	float v = 0.5f - v3_dot(*normal, scene->light) * 0.5f;

	// cheap gamma adjust
	//v = v * v;

	if (style != Draw_Bluenoise)
	{
		int vi = (int)(32.99f * v);
		if (vi > 32)
			vi = 32;
		else if (vi < 0)
			vi = 0;

		// fill
		{
			const uint8_t* pattern = (const uint8_t*)&patterns[vi];
			fillTriangle(bitmap, rowstride, p1, p2, p3, pattern);
		}
		// wire
		if (wire)
		{
			//const uint8_t* color = patterns[0]; // 32: white, 0: black
			vi -= 16;
			if (vi < 0) vi = 0;
			const uint8_t* pattern = (const uint8_t*)&patterns[vi];
			drawLine(bitmap, rowstride, p1, p2, 1, pattern);
			drawLine(bitmap, rowstride, p2, p3, 1, pattern);
			drawLine(bitmap, rowstride, p3, p1, 1, pattern);
		}
	}
	else
	{
		// draw into byte buffer for blue noise thresholding
		int col = (int)(v * 255.0f);
		if (col < 0) col = 0;
		if (col > 255) col = 255;

		draw_triangle(p1, p2, p3, col);

		if (wire)
		{
			col -= 128;
			if (col < 0) col = 0;
			int p1x = (int)p1->x, p1y = (int)p1->y;
			int p2x = (int)p2->x, p2y = (int)p2->y;
			int p3x = (int)p3->x, p3y = (int)p3->y;
			draw_line(g_screen_buffer, SCREEN_X, SCREEN_Y, p1x, p1y, p2x, p2y, col);
			draw_line(g_screen_buffer, SCREEN_X, SCREEN_Y, p2x, p2y, p3x, p3y, col);
			draw_line(g_screen_buffer, SCREEN_X, SCREEN_Y, p3x, p3y, p1x, p1y, col);
		}
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
		drawShapeFace(scene, buffer, rowstride, &scene->tmp_points[idx0], &scene->tmp_points[idx1], &scene->tmp_points[idx2], &scene->tmp_face_normals[fi], style, wire);
	}
}
