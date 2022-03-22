#pragma once

#ifndef VEC_H
#define VEC_H

struct vec_t {
	double x, y;
};
typedef struct vec_t vec_t;

static const vec_t VEC_ZERO;

vec_t vec_from_points(const vec_t* from, const vec_t* to);
double vec_dot(const vec_t* lhs, const vec_t* rhs);
double vec_dist_sq(const vec_t* a, const vec_t* b);
double vec_dist(const vec_t* a, const vec_t* b);
double vec_length_sq(const vec_t* vec);
double vec_length(const vec_t* vec);
void vec_set_length(vec_t* vec, double length);
void vec_normalize(vec_t* vec);

#endif // VEC_H