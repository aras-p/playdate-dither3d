// SPDX-License-Identifier: Unlicense

#include "globals.h"
#include "platform.h"
#include "mathlib.h"
#include "util/pixel_ops.h"
#include "render.h"

#include "fx_meshes.h"

#include <stdlib.h>

#include "gen_scene_data.h"

static Scene s_scene;

static bool s_draw_wire = true;
static enum DrawStyle s_draw_style = Draw_Pattern;

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
	if (G.buttons_cur & kPlatButtonLeft)
	{
		//planeCount -= 1;
		//if (planeCount < 1)
		//	planeCount = 1;
	}
	if (G.buttons_cur & kPlatButtonRight)
	{
		//planeCount += 1;
		//if (planeCount > MAX_PLANES)
		//	planeCount = MAX_PLANES;
	}
	if (G.buttons_pressed & kPlatButtonUp)
	{
		s_draw_style--;
		if (s_draw_style < 0)
			s_draw_style = Draw_Count - 1;
	}
	if (G.buttons_pressed & kPlatButtonDown)
	{
		s_draw_style++;
		if (s_draw_style >= Draw_Count)
			s_draw_style = 0;
	}
	if (G.buttons_pressed & kPlatButtonA)
		s_draw_wire = !s_draw_wire;

	float cangle = G.crank_angle_rad;
	float cs = cosf(cangle);
	float ss = sinf(cangle);
	scene_setCamera(&s_scene, (float3) { cs * 8.0f, 3.0f, ss * 8.0f }, (float3) { 0, 0, 0 }, 1.0f, (float3) { 0, -1, 0 });

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
	if (s_draw_style == Draw_Bluenoise)
		clear_screen_buffers();
	else
		plat_gfx_clear(kSolidColorWhite);

	for (int i = 0; i < SCENE_OBJECT_COUNT; ++i)
	{
		int idx = s_object_order[i];
		scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, g_meshes[idx].mesh, &g_meshes[idx].tr, s_draw_style, s_draw_wire);
	}

	if (s_draw_style == Draw_Bluenoise)
		draw_dithered_screen(G.framebuffer, 0);

	G.statval1 = SCENE_OBJECT_COUNT;
	G.statval2 = s_draw_style;
}

void fx_meshes_init()
{
	scene_init(&s_scene);
	scene_setLight(&s_scene, (float3) { 0.3f, 1.0f, 0.3f });
}
