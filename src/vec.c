#include "vec.h"

#include <math.h>

const vec_t VEC_ZERO = { .x = 0, .y = 0 };
const vec_t VEC_UNIT_X = { .x = 1, .y = 0 };
const vec_t VEC_UNIT_Y = { .x = 0, .y = 1 };

vec_t vec_from_points(const vec_t* from, const vec_t* to) {
	const double dx = to->x - from->x;
	const double dy = to->y - from->y;

	return (vec_t){ .x = dx, .y = dy };
}

double vec_dot(const vec_t* lhs, const vec_t* rhs) {
	return lhs->x * rhs->x + lhs->y * rhs->y;
}

double vec_dist_sq(const vec_t* a, const vec_t* b) {
	const vec_t vec = vec_from_points(a, b);
	return vec_length_sq(&vec);
}

double vec_dist(const vec_t* a, const vec_t* b) {
	return sqrt(vec_dist_sq(a, b));
}

double vec_length_sq(const vec_t* vec) {
	return vec_dot(vec, vec);
}

double vec_length(const vec_t* vec) {
	return sqrt(vec_length_sq(vec));
}

void vec_set_length(vec_t* vec, double length) {
	const double actual_length = vec_length(vec);

	if (fabs(actual_length) > 1e-10) {
		vec->x *= length / actual_length;
		vec->y *= length / actual_length;
	}
}

void vec_normalize(vec_t* vec) {
	vec_set_length(vec, 1.0);
}