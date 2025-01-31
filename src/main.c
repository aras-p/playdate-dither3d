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
	G.ending = false;

	init_pixel_ops();
	fx_plasma_init();
	fx_raytrace_init();
	fx_starfield_init();
	fx_prettyhip_init();
}

static int s_beat_frame_done = -1;
#define TIME_SCRUB_SECONDS (5.0f)

static int track_current_time()
{
	G.frame_count++;
	G.prev_time = G.time;
	G.time = plat_time_get() / TIME_UNIT_LENGTH_SECONDS;

	if (G.prev_time > G.time)
		G.prev_time = G.time;

	// "beat" is if during this frame the tick would change (except when music is done, no beats then)
	int beat_at_end_of_frame = (int)(G.time + TIME_LEN_30FPSFRAME);
	G.beat = (G.ending || (s_beat_frame_done >= beat_at_end_of_frame)) ? false : true;
	return beat_at_end_of_frame;
}

typedef struct DemoEffect {
	float start_time;
	float end_time;
	fx_update_function update;
	float ending_alpha;
} DemoEffect;

static DemoEffect s_effects[] = {
	{0, 32, fx_starfield_update},
	{32, 64, fx_prettyhip_update},
	{64, 96, fx_plasma_update},
	{96, 240, fx_raymarch_update},
	{240, 304, fx_raytrace_update},
};
#define DEMO_EFFECT_COUNT (sizeof(s_effects)/sizeof(s_effects[0]))

static DemoEffect s_ending_effects[] = {
	{0, 32, fx_starfield_update, 0.5f},
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

static int s_cur_effect = DEMO_ENDING_EFFECT_COUNT - 1;


static void update_effect()
{
	if (!G.ending)
	{
		// regular demo part: timeline of effects
		float t = G.time;
		for (int i = 0; i < DEMO_EFFECT_COUNT; ++i)
		{
			const DemoEffect* fx = &s_effects[i];
			if (t >= fx->start_time && t < fx->end_time)
			{
				float a = invlerp(fx->start_time, fx->end_time, t);
				fx->update(fx->start_time, fx->end_time, a);
				break;
			}
		}
	}
	else
	{
		// after demo finish: interactive part
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
}

void app_update()
{
	// track inputs and time
	PlatButtons btCur, btPushed, btRel;
	plat_input_get_buttons(&btCur, &btPushed, &btRel);
	G.buttons_cur = btCur;
	G.buttons_pressed = btPushed;
	G.crank_angle_rad = plat_input_get_crank_angle_rad();

	int beat_at_end_of_frame = track_current_time();

	G.framebuffer = plat_gfx_get_frame();
	G.framebuffer_stride = SCREEN_STRIDE_BYTES;

	// update the effect
	update_effect();

	s_beat_frame_done = beat_at_end_of_frame;

	// draw FPS, time
#if SHOW_STATS
	plat_gfx_draw_stats(G.time);
#endif

	// tell OS that we've updated the whole screen
	plat_gfx_mark_updated_rows(0, SCREEN_Y-1);
}
