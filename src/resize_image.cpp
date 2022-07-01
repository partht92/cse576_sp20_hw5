#include <cmath>
#include <iostream>

#include "image.h"

using namespace std;

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the nearest neighbor to pixel (x,y,c)
float Image::pixel_nearest(float x, float y, int c) const {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)

  int nearest_x = round(x);
  int nearest_y = round(y);

  return clamped_pixel(nearest_x, nearest_y, c);
}

// HW1 #1
// float x,y: inexact coordinates
// int c: channel
// returns the bilinearly interpolated pixel (x,y,c)
float Image::pixel_bilinear(float x, float y, int c) const {
  // Since you are inside class Image you can
  // use the member function pixel(a,b,c)

  int x_low = floor(x);
  int x_high = ceil(x);
  int y_low = floor(y);
  int y_high = ceil(y);

  float A1 = (x - x_low) * (y - y_low);
  float A2 = (x - x_low) * (y_high - y);
  float A3 = (x_high - x) * (y - y_low);
  float A4 = (x_high - x) * (y_high - y);

  return A4 * clamped_pixel(x_low, y_low, c) +
         A3 * clamped_pixel(x_low, y_high, c) +
         A2 * clamped_pixel(x_high, y_low, c) +
         A1 * clamped_pixel(x_high, y_high, c);
}

// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image nearest_resize(const Image& im, int w, int h) {
  Image ret(w, h, im.c);

  float x_scale = (1.0 * w) / im.w;
  float y_scale = (1.0 * h) / im.h;

  auto original_x = [&](int x) { return ((x + 0.5) / x_scale) - 0.5; };
  auto original_y = [&](int y) { return ((y + 0.5) / y_scale) - 0.5; };

  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      // (x’,y’) = (-1/2, -1/2) + scale*(x+1/2,y+1/2)
      for (int c = 0; c < im.c; c++) {
        ret(i, j, c) = im.pixel_nearest(original_x(i), original_y(j), c);
      }
    }
  }

  return ret;
}

// HW1 #1
// int w,h: size of new image
// const Image& im: input image
// return new Image of size (w,h,im.c)
Image bilinear_resize(const Image& im, int w, int h) {
  Image ret(w, h, im.c);

  float x_scale = (1.0 * w) / im.w;
  float y_scale = (1.0 * h) / im.h;

  auto original_x = [&](int x) { return ((x + 0.5) / x_scale) - 0.5; };
  auto original_y = [&](int y) { return ((y + 0.5) / y_scale) - 0.5; };

  for (int i = 0; i < w; i++) {
    for (int j = 0; j < h; j++) {
      // (x’,y’) = (-1/2, -1/2) + scale*(x+1/2,y+1/2)
      for (int c = 0; c < im.c; c++) {
        ret(i, j, c) = im.pixel_bilinear(original_x(i), original_y(j), c);
      }
    }
  }

  return ret;
}
