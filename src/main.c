// SPDX-License-Identifier: Unlicense

#include "platform.h"

#include "effects/fx.h"
#include "globals.h"
#include "mathlib.h"
#include "util/pixel_ops.h"

//#define SHOW_STATS 1

void app_initialize()
{
	G.rng = 1;
	G.frame_count = 0;
	G.time = G.prev_time = -1.0f;

	init_pixel_ops();
	fx_meshes_init();
	fx_plasma_init();
	fx_raytrace_init();
	fx_starfield_init();
	fx_prettyhip_init();
}

#define TIME_SCRUB_SECONDS (5.0f)

static void track_current_time()
{
	G.frame_count++;
	G.prev_time = G.time;
	G.time = plat_time_get() / TIME_UNIT_LENGTH_SECONDS;

	if (G.prev_time > G.time)
		G.prev_time = G.time;
}

typedef struct DemoEffect {
	float start_time;
	float end_time;
	fx_update_function update;
	float ending_alpha;
} DemoEffect;

static DemoEffect s_ending_effects[] = {
	{0, 32, fx_meshes_update, 0.5f},
	{32, 64, fx_prettyhip_update, 0.5f},
	{64, 80, fx_plasma_update, 0.4f}, // w/ twisty cube
	{80, 96, fx_plasma_update, 0.6f}, // w/ twisted ring
	{96, 240, fx_raymarch_update, 0.2f}, // xor towers
	{96, 240, fx_raymarch_update, 0.3f}, // sponge
	{96, 240, fx_raymarch_update, 0.4f}, // puls
	{96, 240, fx_raymarch_update, 0.8f}, // 4x fx rotating
	{240, 304, fx_raytrace_update, 0.5f},
};
#define DEMO_ENDING_EFFECT_COUNT (sizeof(s_ending_effects)/sizeof(s_ending_effects[0]))

static int s_cur_effect = 0;


static void update_effect()
{
	// A/B switch between effects
	if (G.buttons_pressed & kPlatButtonB) {
		clear_screen_buffers();
		s_cur_effect--;
		if (s_cur_effect < 0 || s_cur_effect >= DEMO_ENDING_EFFECT_COUNT)
			s_cur_effect = DEMO_ENDING_EFFECT_COUNT - 1;
	}
	if (G.buttons_pressed & kPlatButtonA) {
		clear_screen_buffers();
		s_cur_effect++;
		if (s_cur_effect < 0 || s_cur_effect >= DEMO_ENDING_EFFECT_COUNT)
			s_cur_effect = 0;
	}
	const DemoEffect* fx = &s_ending_effects[s_cur_effect];
	fx->update(fx->start_time, fx->end_time, fx->ending_alpha);
}

void app_update()
{
	// track inputs and time
	PlatButtons btCur, btPushed, btRel;
	plat_input_get_buttons(&btCur, &btPushed, &btRel);
	G.buttons_cur = btCur;
	G.buttons_pressed = btPushed;
	G.crank_angle_rad = plat_input_get_crank_angle_rad();

	G.framebuffer = plat_gfx_get_frame();
	G.framebuffer_stride = SCREEN_STRIDE_BYTES;

	// update the effect
	update_effect();

	// draw FPS, time
#if SHOW_STATS
	plat_gfx_draw_stats(G.time);
#endif

	// tell OS that we've updated the whole screen
	plat_gfx_mark_updated_rows(0, SCREEN_Y-1);
}
