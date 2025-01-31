// SPDX-License-Identifier: Unlicense

#include "../globals.h"

#include "../platform.h"
#include "fx.h"
#include "../mathlib.h"
#include "../util/pixel_ops.h"
#include "../mini3d/render.h"

#include <stdlib.h>

float3 g_mesh_Cube_vb[] = { // 28 verts
  {1.000000f, -1.000000f, 1.000000f},
  {-1.000000f, -1.000000f, -1.000000f},
  {-0.754643f, 0.700373f, 1.000000f},
  {-0.654210f, -0.654210f, 2.619101f},
  {-0.480073f, 0.111070f, 2.879520f},
  {0.480073f, 0.111070f, 2.879520f},
  {0.654210f, -0.654210f, 2.619101f},
  {-1.000000f, -1.000000f, 1.000000f},
  {0.754643f, 0.408254f, -1.529500f},
  {1.000000f, -1.000000f, -1.000000f},
  {0.754643f, 0.700373f, 1.000000f},
  {-0.754643f, 0.408254f, -1.529500f},
  {-1.283430f, -0.567080f, -1.000000f},
  {-1.283430f, -0.567080f, 1.000000f},
  {1.283430f, 0.050230f, -1.000000f},
  {1.283430f, -0.567080f, -1.000000f},
  {-1.283430f, 0.050230f, 1.000000f},
  {1.283430f, 0.050230f, 1.000000f},
  {1.283430f, -0.567080f, 1.000000f},
  {-1.283430f, 0.050230f, -1.000000f},
  {-3.291110f, -0.774914f, -1.192131f},
  {-3.291110f, -0.774914f, -0.118306f},
  {3.291110f, -0.447044f, -1.192131f},
  {3.291110f, -0.774914f, -1.192131f},
  {-3.291110f, -0.447044f, -0.118306f},
  {3.291110f, -0.447044f, -0.118306f},
  {3.291110f, -0.774914f, -0.118306f},
  {-3.291110f, -0.447044f, -1.192131f},
};
uint16_t g_mesh_Cube_ib[] = { // 52 tris
  8, 11, 2,
  8, 2, 10,
  1, 9, 0,
  1, 0, 7,
  1, 11, 8,
  1, 8, 9,
  6, 5, 4,
  6, 4, 3,
  2, 16, 13,
  2, 13, 4,
  7, 0, 6,
  7, 6, 3,
  10, 2, 4,
  10, 4, 5,
  2, 11, 19,
  2, 19, 16,
  9, 8, 14,
  9, 14, 15,
  4, 13, 7,
  4, 7, 3,
  14, 17, 25,
  14, 25, 22,
  18, 15, 23,
  18, 23, 26,
  11, 1, 12,
  11, 12, 19,
  0, 18, 5,
  0, 5, 6,
  1, 7, 13,
  1, 13, 12,
  0, 9, 15,
  0, 15, 18,
  8, 10, 17,
  8, 17, 14,
  18, 17, 10,
  18, 10, 5,
  21, 24, 27,
  21, 27, 20,
  23, 22, 25,
  23, 25, 26,
  17, 18, 26,
  17, 26, 25,
  15, 14, 22,
  15, 22, 23,
  16, 19, 27,
  16, 27, 24,
  19, 12, 20,
  19, 20, 27,
  13, 16, 24,
  13, 24, 21,
  12, 13, 21,
  12, 21, 20,
};

static Scene s_scene;
static Mesh s_shape_plane;

#define MAX_PLANES 500
static int planeCount = 100;
static xform planeXforms[MAX_PLANES];
static float planeDistances[MAX_PLANES];
static int planeOrder[MAX_PLANES];

static int CompareZ(const void* a, const void* b)
{
	int ia = *(const int*)a;
	int ib = *(const int*)b;
	float za = planeDistances[ia];
	float zb = planeDistances[ib];
	if (za < zb)
		return +1;
	if (za > zb)
		return -1;
	return 0;
}


void fx_meshes_update(float start_time, float end_time, float alpha)
{
	if (G.buttons_cur & kPlatButtonLeft)
	{
		planeCount -= 1;
		if (planeCount < 1)
			planeCount = 1;
	}
	if (G.buttons_cur & kPlatButtonRight)
	{
		planeCount += 1;
		if (planeCount > MAX_PLANES)
			planeCount = MAX_PLANES;
	}

	uint32_t rng = 1;

	float cangle = G.crank_angle_rad;
	float cs = cosf(cangle);
	float ss = sinf(cangle);
	scene_setCamera(&s_scene, (float3) { cs * 8.0f, 3.0f, ss * 8.0f }, (float3) { 0, 0, 0 }, 1.0f, (float3) { 0, -1, 0 });

	// position and sort planes
	for (int i = 0; i < planeCount; ++i)
	{
		float px = ((XorShift32(&rng) & 63) - 31.5f);
		float py = ((XorShift32(&rng) & 15) - 7.5f);
		float pz = ((XorShift32(&rng) & 63) - 31.5f);
		float rot = (XorShift32(&rng) % 360) * 1.0f;
		float rx = ((XorShift32(&rng) & 63) - 31.5f);
		float ry = 60;
		float rz = ((XorShift32(&rng) & 63) - 31.5f);
		if (i == 0) {
			px = py = pz = 0.0f;
			rot = 0;
			rx = 0;
			ry = 1;
			rz = 0;
		}
		planeXforms[i] = xform_make_axis_angle(rot, (float3) { rx, ry, rz });
		planeXforms[i].x = px;
		planeXforms[i].y = py;
		planeXforms[i].z = pz;
		//planeXforms[i] = mtx_make_translate(px, py, pz);
		float3 center = xform_transform_pt(&s_scene.camera, xform_transform_pt(&planeXforms[i], (float3) {0,0,0}));
		planeDistances[i] = center.z;
		planeOrder[i] = i;
	}
	qsort(planeOrder, planeCount, sizeof(planeOrder[0]), CompareZ);

	// draw
	plat_gfx_clear(kSolidColorWhite);
	for (int i = 0; i < planeCount; ++i)
	{
		int idx = planeOrder[i];
		scene_drawMesh(&s_scene, G.framebuffer, G.framebuffer_stride, &s_shape_plane, &planeXforms[idx]);
	}
}

void fx_meshes_init()
{
	scene_init(&s_scene);
	scene_setLight(&s_scene, (float3) { 0.3f, 1.0f, 0.3f });

	s_shape_plane.vertices = g_mesh_Cube_vb;
	s_shape_plane.vertex_count = sizeof(g_mesh_Cube_vb) / sizeof(g_mesh_Cube_vb[0]);
	s_shape_plane.tris = g_mesh_Cube_ib;
	s_shape_plane.tri_count = sizeof(g_mesh_Cube_ib) / sizeof(g_mesh_Cube_ib[0]) / 3;
}
