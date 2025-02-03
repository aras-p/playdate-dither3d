# Surface-Stable Fractal Dithering on Playdate

Rune Skovbo Johansen has a really sweet [**Surface-Stable Fractal Dithering**](https://github.com/runevision/Dither3D)
technique, where the dither dots "stick" to 3D surfaces, yet the dot density adapts to the view distance and zoom
level.

Some people have asked whether this would be a good technique for [Playdate](https://play.date/), given that the screen
is one-bit color. And so I had to try it out!

**My impression: not practical, really**. Playdate hardware is like a PC from 1995 - no GPU at all, one fairly simple CPU
core. As such, it can do fairly simple 3D rendering (well, you need to write the whole rasterizer on the CPU),
but can barely do more than a handful of math operations per rasterized pixel. Rasterizing with screen-space fixed
Bayer or Blue Noise dither patterns is the way to go due to their simplicity.

![Screenshot](/img/250203b-pd-interp_x1.png?raw=true "Screenshot")

However! It was fairly fun hacking on this. Amusingly enough, I have never written a perspective-correct interpolating
triangle rasterizer, etc. So this was a good learning experience. For reference, a non-perspective-correct vs. correct
texture coordinate interpolation:
|Not correct |Correct |
|-|-|
| ![Screenshot](/img/250131b-checker-side.png?raw=true "Screenshot") | ![Screenshot](/img/250131c-checker-side-persp.png?raw=true "Screenshot") |

Current status:
- The really simple scene in the first image runs at **13 frames per second** on a Playdate. Initially I had it running at 2FPS.
  Currently it has some simplifications and behavior changes compared to full Rune's technique:
  - Dot spacing is fixed (6.0),
  - Dither pattern is 4x4 Bayer one, but with texture XY resolution reduced twice (i.e. 3D texture is 32x32x16).
  - No anisotropic derivatives handling nor contrast tweaking based on that.
- Rasterizer is a half-space based one, processing 2x2 pixel blocks in one iteration, and doing perspective
  correct UV interpolation every 2 pixels horizontally.
- Very likely that a traditional "scanline" triangle rasterizer would be a better fit for the
  hardware. I should try that out someday.
- It is entirely possible that there's a ton of low hanging fruit w.r.t. optimizations that I have
  overlooked so far.


### License

Everything I wrote myself is Unlicense / Public Domain. However some 3rd party libraries are used too:
- No direct code used, but the whole [**Dither3D technique**](https://github.com/runevision/Dither3D) is MPL-2.0 license.
- (only on PC): [**Sokol**](https://github.com/floooh/sokol): sokol_app, sokol_gfx, sokol_time. zlib/libpng license.
