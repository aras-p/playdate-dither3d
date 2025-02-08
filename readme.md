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

![Screenshot](/img/250208a-pd-dither-scanline.png?raw=true "Screenshot")

However! It was fairly fun hacking on this. Amusingly enough, I have never written a perspective-correct interpolating
triangle rasterizer, etc. So this was a good learning experience. For reference, a non-perspective-correct vs. correct
texture coordinate interpolation:
|Not correct |Correct |
|-|-|
| ![Screenshot](/img/250131b-checker-side.png?raw=true "Screenshot") | ![Screenshot](/img/250131c-checker-side-persp.png?raw=true "Screenshot") |

Current status:
- The really simple scene in the first image runs at **20 frames per second** on a Playdate. Initially I had it running at 1.5FPS.
  Currently it has some simplifications and behavior changes compared to full Rune's technique:
  - Dot spacing is fixed (6.0), contrast is fixed (1.0), size variability too (0.0).
  - Dither pattern is 4x4 Bayer one, but with texture XY resolution reduced twice (i.e. 3D texture is 32x32x16).
  - No anisotropic derivatives handling nor contrast tweaking based on that.
- I have two rasterizer approaches:
  - Traditional "scanline" one, very much like in Chris Hecker's 1995/1996 article series. Perspective correct
    UV interpolation is done every 8 pixels. The dither pattern spacing is calculated only at triangle vertices
    and interpolated across. Runs at 48ms/frame, code in `draw_triangle_dither3d_scanline`.
  - Halfspace/barycentric one, processing 2x2 pixel blocks in one iteration, and doing perspective
    correct UV interpolation every 2 pixels horizontally. Runs at 63ms/frame, code in `draw_triangle_dither3d_halfspace`.
- It is entirely possible that there's a ton of low hanging fruit w.r.t. optimizations that I have
  overlooked so far.

Most of the code layout and structure is taken from [Everybody Wants to Crank the World](https://github.com/aras-p/demo-pd-cranktheworld)
demo I made in 2024. As such, it compiles on Playdate, and on PC as well (via Sokol).

### License

Everything I wrote myself is Unlicense / Public Domain. However some 3rd party libraries are used too:
- No direct code used, but the whole [**Dither3D technique**](https://github.com/runevision/Dither3D) is MPL-2.0 license.
- Scanline rasterizer is based on Chris Hecker's ["Perspective Texture Mapping" series](https://chrishecker.com/Miscellaneous_Technical_Articles) (1995-1996),
  specifically `SUBAFXFL.CPP`.
- Halfspace rasterizer is based on Fabien Giesen's "[Simple watertight triangle rasterizer](https://gist.github.com/rygorous/9b793cd21d876da928bf4c7f3e625908)"
  and [blog posts](https://fgiesen.wordpress.com/2013/02/10/optimizing-the-basic-rasterizer/).
- Parts were based on `mini3d` from Playdate SDK examples (Created by Dave Hayden on 10/20/15, Copyright © 2015 Panic, Inc.)
  and on [mini3d-plus](https://github.com/nstbayless/mini3d-plus) (MIT license), but by now I have replaced them.
- (only on PC): [Sokol](https://github.com/floooh/sokol): sokol_app, sokol_gfx, sokol_time. zlib/libpng license.
