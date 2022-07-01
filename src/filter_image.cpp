#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "image.h"

#define M_PI 3.14159265358979323846

// HW1 #2.1
// Image& im: image to L1-normalize
void l1_normalize(Image& im) {
  float sum = 0;
  for (int i = 0; i < im.w * im.h * im.c; i++) {
    sum += im.data[i];
  }
  for (int i = 0; i < im.w * im.h * im.c; i++) {
    im.data[i] /= sum;
  }
}

// HW1 #2.1
// int w: size of filter
// returns the filter Image of size WxW
Image make_box_filter(int w) {
  assert(w % 2);  // w needs to be odd
  Image box_filter(w, w, 1);

  for (int i = 0; i < w * w; i++) {
    box_filter.data[i] = 1.0 / (w * w);
  }
  l1_normalize(box_filter);

  return box_filter;
}

// HW1 #2.2
// const Image&im: input image
// const Image& filter: filter to convolve with
// bool preserve: whether to preserve number of channels
// returns the convolved image
Image convolve_image(const Image& im, const Image& filter, bool preserve) {
  assert(filter.c == 1);

  int channels = preserve ? im.c : 1;
  Image ret(im.w, im.h, channels);
  // This is the case when we need to use the function clamped_pixel(x,y,c).
  // Otherwise you'll have to manually check whether the filter goes out of
  // bounds

  // TODO: Make sure you set the sizes of ret properly. Use ret=Image(w,h,c) to
  // reset ret
  // TODO: Do the convolution operator
  auto convolve = [&](int x, int y, int c) {
    float value = 0;
    for (int i = 0; i < filter.w; i++) {
      for (int j = 0; j < filter.h; j++) {
        value +=
            im.clamped_pixel(x + i - filter.w / 2, y + j - filter.h / 2, c) *
            filter(i, j, 0);
      }
    }
    return value;
  };

  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      float agg_value = 0;
      for (int c = 0; c < channels; c++) {
        if (preserve) {
          ret.set_pixel(i, j, c, convolve(i, j, c));
        } else {
          agg_value += convolve(i, j, c);
        }
      }
      if (!preserve) ret.set_pixel(i, j, 0, agg_value);
    }
  }

  return ret;
}

// HW1 #2.3
// returns basic 3x3 high-pass filter
Image make_highpass_filter() {
  Image high_pass_filter(3, 3, 1);
  high_pass_filter.set_pixel(0, 0, 0, 0);
  high_pass_filter.set_pixel(1, 0, 0, -1);
  high_pass_filter.set_pixel(2, 0, 0, 0);
  high_pass_filter.set_pixel(0, 1, 0, -1);
  high_pass_filter.set_pixel(1, 1, 0, 4);
  high_pass_filter.set_pixel(2, 1, 0, -1);
  high_pass_filter.set_pixel(0, 2, 0, 0);
  high_pass_filter.set_pixel(1, 2, 0, -1);
  high_pass_filter.set_pixel(2, 2, 0, 0);
  return high_pass_filter;
}

// HW1 #2.3
// returns basic 3x3 sharpen filter
Image make_sharpen_filter() {
  Image sharpen_filter(3, 3, 1);
  sharpen_filter.set_pixel(0, 0, 0, 0);
  sharpen_filter.set_pixel(1, 0, 0, -1);
  sharpen_filter.set_pixel(2, 0, 0, 0);
  sharpen_filter.set_pixel(0, 1, 0, -1);
  sharpen_filter.set_pixel(1, 1, 0, 5);
  sharpen_filter.set_pixel(2, 1, 0, -1);
  sharpen_filter.set_pixel(0, 2, 0, 0);
  sharpen_filter.set_pixel(1, 2, 0, -1);
  sharpen_filter.set_pixel(2, 2, 0, 0);
  return sharpen_filter;
}

// HW1 #2.3
// returns basic 3x3 emboss filter
Image make_emboss_filter() {
  Image emboss_filter(3, 3, 1);

  emboss_filter.set_pixel(0, 0, 0, -2);
  emboss_filter.set_pixel(1, 0, 0, -1);
  emboss_filter.set_pixel(2, 0, 0, 0);
  emboss_filter.set_pixel(0, 1, 0, -1);
  emboss_filter.set_pixel(1, 1, 0, 1);
  emboss_filter.set_pixel(2, 1, 0, 1);
  emboss_filter.set_pixel(0, 2, 0, 0);
  emboss_filter.set_pixel(1, 2, 0, 1);
  emboss_filter.set_pixel(2, 2, 0, 2);

  return emboss_filter;
}

// HW1 #2.4
// float sigma: sigma for the gaussian filter
// returns basic gaussian filter
Image make_gaussian_filter(float sigma) {
  int six_sigma = ceil(6 * sigma);
  int filter_dimension = six_sigma % 2 == 0 ? six_sigma + 1 : six_sigma;

  Image filter(filter_dimension, filter_dimension, 1);
  for (int i = 0; i < filter_dimension; i++) {
    for (int j = 0; j < filter_dimension; j++) {
      int x_dist = i - filter_dimension / 2;
      int y_dist = j - filter_dimension / 2;
      filter.set_pixel(
          i, j, 0,
          exp(-1 * (x_dist * x_dist + y_dist * y_dist) / (2 * sigma * sigma)) /
              (2 * M_PI * sigma * sigma));
    }
  }
  l1_normalize(filter);
  return filter;
}

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their sum
Image add_image(const Image& a, const Image& b) {
  assert(a.w == b.w && a.h == b.h &&
         a.c == b.c);  // assure images are the same size

  Image added_image(a.w, a.h, a.c);
  for (int i = 0; i < a.w * a.c * a.h; i++) {
    added_image.data[i] = a.data[i] + b.data[i];
  }

  return added_image;
}

// HW1 #3
// const Image& a: input image
// const Image& b: input image
// returns their difference res=a-b
Image sub_image(const Image& a, const Image& b) {
  assert(a.w == b.w && a.h == b.h &&
         a.c == b.c);  // assure images are the same size

  Image subtracted_image(a.w, a.h, a.c);
  for (int i = 0; i < a.w * a.c * a.h; i++) {
    subtracted_image.data[i] = a.data[i] - b.data[i];
  }

  return subtracted_image;
}

// HW1 #4.1
// returns basic GX filter
Image make_gx_filter() {
  Image gx_filter(3, 3, 1);

  gx_filter.set_pixel(0, 0, 0, -1);
  gx_filter.set_pixel(1, 0, 0, 0);
  gx_filter.set_pixel(2, 0, 0, 1);
  gx_filter.set_pixel(0, 1, 0, -2);
  gx_filter.set_pixel(1, 1, 0, 0);
  gx_filter.set_pixel(2, 1, 0, 2);
  gx_filter.set_pixel(0, 2, 0, -1);
  gx_filter.set_pixel(1, 2, 0, 0);
  gx_filter.set_pixel(2, 2, 0, 1);

  return gx_filter;
}

// HW1 #4.1
// returns basic GY filter
Image make_gy_filter() {
  Image gy_filter(3, 3, 1);

  gy_filter.set_pixel(0, 0, 0, -1);
  gy_filter.set_pixel(1, 0, 0, -2);
  gy_filter.set_pixel(2, 0, 0, -1);
  gy_filter.set_pixel(0, 1, 0, 0);
  gy_filter.set_pixel(1, 1, 0, 0);
  gy_filter.set_pixel(2, 1, 0, 0);
  gy_filter.set_pixel(0, 2, 0, 1);
  gy_filter.set_pixel(1, 2, 0, 2);
  gy_filter.set_pixel(2, 2, 0, 1);

  return gy_filter;
}

// HW1 #4.2
// Image& im: input image
void feature_normalize(Image& im) {
  assert(im.w * im.h);  // assure we have non-empty image
  for (int c = 0; c < im.c; c++) {
    float max = im(0, 0, c);
    float min = im(0, 0, c);
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        min = std::min(min, im(i, j, c));
        max = std::max(max, im(i, j, c));
      }
    }

    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        im.set_pixel(i, j, c,
                     (max - min == 0) ? 0 : (im(i, j, c) - min) / (max - min));
      }
    }
  }
}

// Normalizes features across all channels
void feature_normalize_total(Image& im) {
  assert(im.w * im.h * im.c);  // assure we have non-empty image

  int nc = im.c;
  im.c = 1;
  im.w *= nc;

  feature_normalize(im);

  im.w /= nc;
  im.c = nc;
}

// HW1 #4.3
// Image& im: input image
// return a pair of images of the same size
pair<Image, Image> sobel_image(const Image& im) {
  Image gx_image = convolve_image(im, make_gx_filter(), /*preserve=*/false);
  Image gy_image = convolve_image(im, make_gy_filter(), /*preserve=*/false);

  Image magnitude(im.w, im.h, 1);
  Image direction(im.w, im.h, 1);
  for (int i = 0; i < im.w; i++) {
    for (int j = 0; j < im.h; j++) {
      float gx = gx_image(i, j, 0);
      float gy = gy_image(i, j, 0);
      magnitude.set_pixel(i, j, 0, std::sqrtf(gx * gx + gy * gy));
      direction.set_pixel(i, j, 0, std::atan2f(gy, gx) / (2 * M_PI) + 0.5);
    }
  }
  feature_normalize(magnitude);

  return {magnitude, direction};
}

// HW1 #4.4
// const Image& im: input image
// returns the colorized Sobel image of the same size
Image colorize_sobel(const Image& im) {
  // TODO: Your code here
  NOT_IMPLEMENTED();

  return im;
}

// HW1 #4.5
// const Image& im: input image
// float sigma1,sigma2: the two sigmas for bilateral filter
// returns the result of applying bilateral filtering to im
Image bilateral_filter(const Image& im, float sigma1, float sigma2) {
  Image bf(im.w, im.h, im.c);

  const auto gaussian = [&](float x_dist, float sigma) {
    return exp(-1 * (x_dist * x_dist) / (2 * sigma * sigma)) /
           (2 * M_PI * sigma * sigma);
  };

  int six_sigma = ceil(6 * sigma1);
  int filter_dimension = six_sigma % 2 == 0 ? six_sigma + 1 : six_sigma;

  const auto apply_filter = [&](int x, int y, int c) {
    float value = 0;
    float normalize_value = 0;
    for (int i = 0; i < filter_dimension; i++) {
      for (int j = 0; j < filter_dimension; j++) {
        int x_prime = x + i - filter_dimension / 2;
        int y_prime = y + j - filter_dimension / 2;

        float weight = gaussian(x - x_prime, sigma1) *
                       gaussian(y - y_prime, sigma1) *
                       gaussian(im.clamped_pixel(x, y, c) -
                                    im.clamped_pixel(x_prime, y_prime, c),
                                sigma2);

        normalize_value += weight;
        value += im.clamped_pixel(x_prime, y_prime, c) * weight;
      }
    }

    return value / normalize_value;
  };

  for (int c = 0; c < im.c; c++) {
    for (int i = 0; i < im.w; i++) {
      for (int j = 0; j < im.h; j++) {
        bf.set_pixel(i, j, c, apply_filter(i, j, c));
      }
    }
  }

  return bf;
}

// HELPER MEMBER FXNS

void Image::feature_normalize(void) { ::feature_normalize(*this); }
void Image::feature_normalize_total(void) { ::feature_normalize_total(*this); }
void Image::l1_normalize(void) { ::l1_normalize(*this); }

Image operator-(const Image& a, const Image& b) { return sub_image(a, b); }
Image operator+(const Image& a, const Image& b) { return add_image(a, b); }
