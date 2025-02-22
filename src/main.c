// SPDX-License-Identifier: Unlicense

#include "platform.h"

#include "fx_meshes.h"
#include "globals.h"
#include "mathlib.h"
#include "util/pixel_ops.h"

#define SHOW_STATS 1

void app_initialize()
{
	init_pixel_ops();
	fx_meshes_init();
}

void app_update()
{
#if SHOW_STATS
	plat_gfx_begin_stats();
#endif

	// track inputs
	PlatButtons btCur, btPushed, btRel;
	plat_input_get_buttons(&btCur, &btPushed, &btRel);
	G.buttons_cur = btCur;
	G.buttons_pressed = btPushed;
	G.crank_angle_rad = plat_input_get_crank_angle_rad();

	G.framebuffer = plat_gfx_get_frame();
	G.framebuffer_stride = SCREEN_STRIDE_BYTES;

	// update the effect
	fx_meshes_update();

	// draw FPS, stats
#if SHOW_STATS
	plat_gfx_draw_stats(G.statval1, G.statval2);
#endif

	// tell OS that we've updated the whole screen
	plat_gfx_mark_updated_rows(0, SCREEN_Y-1);
}
