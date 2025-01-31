// SPDX-License-Identifier: Unlicense

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

#define TIME_UNIT_LENGTH_SECONDS (0.46875f)
#define TIME_LEN_30FPSFRAME (1.0f/(TIME_UNIT_LENGTH_SECONDS*30.0f))

typedef struct Globals
{
	// time: all units in music beats (1.0 = 1 beat = 0.46875s; 4.0 = 1 bar = 1.875s)
	float time; // current global time
	float prev_time;
	int frame_count;

	// screen
	uint8_t* framebuffer;
	int framebuffer_stride;

	// misc
	uint32_t rng;
	int statval1;
	int statval2;

	// input
	uint32_t buttons_cur; // buttons currently pressed, bitmask of kButton*
	uint32_t buttons_pressed; // buttons that were pressed down this frame, bitmask of kButton*
	float crank_angle_rad;
} Globals;

extern Globals G;
