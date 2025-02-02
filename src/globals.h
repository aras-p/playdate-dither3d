// SPDX-License-Identifier: Unlicense

#pragma once

#include <stddef.h>
#include <stdint.h>
#include <stdbool.h>

typedef struct Globals
{
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
