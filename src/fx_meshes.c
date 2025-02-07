// SPDX-License-Identifier: Unlicense

#include "globals.h"
#include "platform.h"
#include "mathlib.h"
#include "util/pixel_ops.h"
#include "render.h"

#include "fx_meshes.h"

#include <stdlib.h>

#include "gen_scene_data.h"

// -------- Test Plane
static float3 g_mesh_TestPlane_vb[] = { // 4 verts
  {-0.5f, -0.5f, 0.0f},
  {-0.5f,  0.5f, 0.0f},
  { 0.5f,  0.5f, 0.0f},
  { 0.5f, -0.5f, 0.0f},
};
static uint16_t g_mesh_TestPlane_ib[] = { // 2 tris
  0, 1, 2,
  0, 2, 3,
};
static float g_mesh_TestPlane_uv[] = { // 6 UV coords
	0, 0,
	0, 2,
	2, 2,
	0, 0,
	2, 2,
	2, 0,
};
static Mesh g_mesh_TestPlane = { 4, g_mesh_TestPlane_vb, 2, g_mesh_TestPlane_ib, g_mesh_TestPlane_uv };

static Scene s_scene;

static bool s_draw_test = false;
static bool s_draw_wire = false;
static enum DrawStyle s_draw_style = Draw_Dither3D_Scanline;
static float s_cam_dist = 8.0f;

#define SCENE_OBJECT_COUNT (sizeof(g_meshes)/sizeof(g_meshes[0]))

static float s_object_distances[SCENE_OBJECT_COUNT];
static int s_object_order[SCENE_OBJECT_COUNT];

static int CompareZ(const void* a, const void* b)
{
	int ia = *(const int*)a;
	int ib = *(const int*)b;
	float za = s_object_distances[ia];
	float zb = s_object_distances[ib];
	if (za < zb)
		return +1;
	if (za > zb)
		return -1;
	return 0;
}

void fx_meshes_update()
{
	if (G.buttons_cur & kPlatButtonUp)
	{
		s_cam_dist *= 0.9f;
		s_cam_dist = max2f(3.0f, s_cam_dist);
	}
	if (G.buttons_cur & kPlatButtonDown)
	{
		s_cam_dist *= 1.1f;
		s_cam_dist = min2f(20.0f, s_cam_dist);
	}
	if (G.buttons_pressed & kPlatButtonLeft)
	{
		s_draw_style--;
		if (s_draw_style < 0 || s_draw_style >= Draw_Count)
			s_draw_style = Draw_Count - 1;
	}
	if (G.buttons_pressed & kPlatButtonRight)
	{
		s_draw_style++;
		if (s_draw_style < 0 || s_draw_style >= Draw_Count)
			s_draw_style = 0;
	}
	if (G.buttons_pressed & kPlatButtonA)
		s_draw_wire = !s_draw_wire;
	if (G.buttons_pressed & kPlatButtonB)
		s_draw_test = !s_draw_test;

	float cangle = G.crank_angle_rad;
	float cs = cosf(cangle);
	float ss = sinf(cangle);
	scene_setCamera(&s_scene, (float3) {
		cs * s_cam_dist, 3.0f, ss * s_cam_dist
	}, (float3) { 0, 0, 0 }, 1.0f, (float3) { 0, -1, 0 });

	// sort the objects
	for (int i = 0; i < SCENE_OBJECT_COUNT; ++i)
	{
		float3 center_local = (float3){ g_meshes[i].tr.x, g_meshes[i].tr.y, g_meshes[i].tr.z };
		float3 center = xform_transform_pt(&s_scene.camera, center_local);
		if (g_meshes[i].mesh == &g_mesh_Ground)
			center.z = 1000.0f;
		s_object_distances[i] = center.z;
		s_object_order[i] = i;
	}
	qsort(s_object_order, SCENE_OBJECT_COUNT, sizeof(s_object_order[0]), CompareZ);

	// draw
	plat_gfx_clear(kSolidColorWhite);

	if (!s_draw_test)
	{
		for (int i = 0; i < SCENE_OBJECT_COUNT; ++i)
		{
			int idx = s_object_order[i];
			scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, g_meshes[idx].mesh, &g_meshes[idx].tr, s_draw_style, s_draw_wire);
		}
	}
	else
	{
		float3 p1, p2, p3, nn;
		nn = f3(0, 0, -1);

		// right side, from the bottom:
		// 20x20 in bottom right, 1px away from screen edge
		p1 = f3(SCREEN_X - 1 - 20, SCREEN_Y - 1 - 20, 1); p2 = f3(SCREEN_X - 1 - 20, SCREEN_Y - 1, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 1, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 1 - 20, SCREEN_Y - 1 - 20, 1); p2 = f3(SCREEN_X - 1, SCREEN_Y - 1, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 1 - 20, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// 30x30, so each checker square maps to 1.5 pixels, should be aliased or non-uniform somehow
		p1 = f3(SCREEN_X - 1 - 30, SCREEN_Y - 2 - 50, 1); p2 = f3(SCREEN_X - 1 - 30, SCREEN_Y - 2 - 20, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 2 - 20, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 1 - 30, SCREEN_Y - 2 - 50, 1); p2 = f3(SCREEN_X - 1, SCREEN_Y - 2 - 20, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 2 - 50, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// 40x40, so each checker square maps to 2
		p1 = f3(SCREEN_X - 1 - 40, SCREEN_Y - 3 - 90, 1); p2 = f3(SCREEN_X - 1 - 40, SCREEN_Y - 3 - 50, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 3 - 50, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 1 - 40, SCREEN_Y - 3 - 90, 1); p2 = f3(SCREEN_X - 1, SCREEN_Y - 3 - 50, 1); p3 = f3(SCREEN_X - 1, SCREEN_Y - 3 - 90, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// left side, from the bottom:

		// 20x20 in bottom left, 1px away from screen edge
		p1 = f3(1, SCREEN_Y-1-20, 1); p2 = f3(1, SCREEN_Y - 1, 1); p3 = f3(1+20, SCREEN_Y - 1, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(1, SCREEN_Y - 1 - 20, 1); p2 = f3(1 + 20, SCREEN_Y - 1, 1); p3 = f3(1 + 20, SCREEN_Y - 1 - 20, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// 20x20 atop of it, 1px gap, rotated 90
		p1 = f3(1 + 20, SCREEN_Y - 2 - 40, 1); p2 = f3(1, SCREEN_Y - 2 - 40, 1); p3 = f3(1, SCREEN_Y - 2 - 20, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(1 + 20, SCREEN_Y - 2 - 40, 1); p2 = f3(1, SCREEN_Y - 2 - 20, 1); p3 = f3(1 + 20, SCREEN_Y - 2 - 20, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// same atop of it, offset by less than half a pixel, should look the same
		p1 = f3(0.6f + 20, SCREEN_Y - 3.3f - 60, 1); p2 = f3(0.6f, SCREEN_Y - 3.3f - 60, 1); p3 = f3(0.6f, SCREEN_Y - 3.3f - 40, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(0.6f + 20, SCREEN_Y - 3.3f - 60, 1); p2 = f3(0.6f, SCREEN_Y - 3.3f - 40, 1); p3 = f3(0.6f + 20, SCREEN_Y - 3.3f - 40, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// clipped by screen edges
		// bottom
		p1 = f3(SCREEN_X / 2 - 10, SCREEN_Y - 30, 1); p2 = f3(SCREEN_X / 2 - 30, SCREEN_Y + 10, 1); p3 = f3(SCREEN_X / 2 + 10, SCREEN_Y + 30, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X / 2 - 10, SCREEN_Y - 30, 1); p2 = f3(SCREEN_X / 2 + 10, SCREEN_Y + 30, 1); p3 = f3(SCREEN_X / 2 + 30, SCREEN_Y - 10, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// top
		p1 = f3(SCREEN_X / 2 + 40 - 10, 0 - 30, 1); p2 = f3(SCREEN_X / 2 + 40 - 30, 0 + 10, 1); p3 = f3(SCREEN_X / 2 + 40 + 10, 0 + 30, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X / 2 + 40 - 10, 0 - 30, 1); p2 = f3(SCREEN_X / 2 + 40 + 10, 0 + 30, 1); p3 = f3(SCREEN_X / 2 + 40 + 30, 0 - 10, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// left
		p1 = f3(0 - 10, SCREEN_Y/2 - 30, 1); p2 = f3(0 - 30, SCREEN_Y/2 + 10, 1); p3 = f3(0 + 10, SCREEN_Y/2 + 30, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(0 - 10, SCREEN_Y/2 - 30, 1); p2 = f3(0 + 10, SCREEN_Y/2 + 30, 1); p3 = f3(0 + 30, SCREEN_Y/2 - 10, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// top right
		p1 = f3(SCREEN_X - 10, 0 - 30, 1); p2 = f3(SCREEN_X - 30, 0 + 10, 1); p3 = f3(SCREEN_X + 10, 0 + 30, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 10, 0 - 30, 1); p2 = f3(SCREEN_X + 10, 0 + 30, 1); p3 = f3(SCREEN_X + 30, 0 - 10, 1);
		drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);

		// quad meshes with perspective
		scene_setCamera(&s_scene, f3(0, 0, -1.5f), f3(0, 0, 0), 1.0f, f3(0, 1, 0));
		xform tr = xform_make_axis_angle(M_PIf * -0.5f, f3(1, 0, 0));
		tr.y = 0.3f;
		scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, &g_mesh_TestPlane, &tr, s_draw_style, s_draw_wire);

		tr = xform_make_axis_angle(M_PIf * 0.5f, f3(1, 0, 0));
		xform rrr = xform_make_axis_angle(M_PIf * 0.3f, f3(0, 1, 0));
		tr = xform_multiply(&tr, &rrr);
		tr.y = -0.3f;
		tr.x = -0.5f;
		scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, &g_mesh_TestPlane, &tr, s_draw_style, s_draw_wire);

		tr = xform_make_axis_angle(M_PIf * -0.8f, f3(0, 0, 1));
		//rrr = xform_make_axis_angle(M_PIf * 0.3f, f3(0, 1, 0));
		//tr = xform_multiply(&tr, &rrr);
		tr.y = 0.3f;
		tr.x = -2.0f;
		tr.z = 2;
		scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, &g_mesh_TestPlane, &tr, s_draw_style, s_draw_wire);

		// top right area, adjoining triangles in various configuration; should be no gaps and no double raster
		const float sc = 2.0f;
		p1 = f3(SCREEN_X - 60 * sc, 20 * sc, 1); p2 = f3(SCREEN_X - 80 * sc, 40 * sc, 1); p3 = f3(SCREEN_X - 60 * sc, 60 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 60 * sc, 20 * sc, 1); p2 = f3(SCREEN_X - 60 * sc, 60 * sc, 1); p3 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); p2 = f3(SCREEN_X - 60 * sc, 60 * sc, 1); p3 = f3(SCREEN_X - 30 * sc, 90 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); p2 = f3(SCREEN_X - 30 * sc, 90 * sc, 1); p3 = f3(SCREEN_X - 10 * sc, 50 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 50 * sc, 10 * sc, 1); p2 = f3(SCREEN_X - 60 * sc, 20 * sc, 1); p3 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 50 * sc, 10 * sc, 1); p2 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); p3 = f3(SCREEN_X - 30 * sc, 10 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 30 * sc, 10 * sc, 1); p2 = f3(SCREEN_X - 40 * sc, 40 * sc, 1); p3 = f3(SCREEN_X - 10 * sc, 50 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 0, s_draw_style, s_draw_wire);
		p1 = f3(SCREEN_X - 30 * sc, 10 * sc, 1); p2 = f3(SCREEN_X - 10 * sc, 50 * sc, 1); p3 = f3(SCREEN_X - 10 * sc, 10 * sc, 1); drawShapeFace(&s_scene, G.framebuffer, G.framebuffer_stride, &p1, &p2, &p3, &nn, &g_mesh_TestPlane, 1, s_draw_style, s_draw_wire);
	}

	G.statval1 = s_draw_style;
	G.statval2 = s_draw_wire ? 1 : 0;
}

void fx_meshes_init()
{
	scene_init(&s_scene);
	scene_setLight(&s_scene, (float3) { -1.0f, -1.0f, 0.9f });
}
