#ifndef COLOR_H
#define COLOR_H

#include "vec3.h"
#include <iostream>

typedef vec3 color;

void write_color(std::ostream& out, const color& pixel_color)
{
  int r = int(255.999 * pixel_color.x);
  int g = int(255.999 * pixel_color.y);
  int b = int(255.999 * pixel_color.z);

  // out << r << ' ' << g << ' ' << b << '\n';
}

#endif