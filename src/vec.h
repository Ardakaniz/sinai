#pragma once

#ifndef VEC_H
#define VEC_H

#define VEC_ZERO (vec_t){ .x = 0, .y = 0 }

struct vec_t {
	double x, y;
};
typedef struct vec_t vec_t;


vec_t vec_from_points(const vec_t* from, const vec_t* to);
double vec_dot(const vec_t* lhs, const vec_t* rhs);
double vec_dist_sq(const vec_t* a, const vec_t* b);
double vec_dist(const vec_t* a, const vec_t* b);
double vec_length_sq(const vec_t* vec);
double vec_length(const vec_t* vec);
void vec_set_length(vec_t* vec, double length);
void vec_normalize(vec_t* vec);

#endif // VEC_H