#ifndef RAY_H
#define RAY_H

#include "vec3.h"

struct ray{
  vec3 ori;
  vec3 dir;

  // constructors
  ray() {}
  ray(const vec3& origin, const vec3& direction) { ori = origin; dir = direction; }

  // functions
  const vec3& origin() const { return ori; }
  const vec3& direction() const { return dir; }

  // P(t) = A + tb
  vec3 at(double t) const { return ori + t * dir; }
};

#endif