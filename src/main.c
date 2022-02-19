#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#include "vec.h"

#define ATOM_RADIUS   0.25 // As a proportion of lattice's size
#define GRID_INTERVAL 0.005
#define GRID_SIZE (size_t)((1.0 / GRID_INTERVAL) * (1.0 / GRID_INTERVAL))

/* STRUCTS */
struct data_t {	
	unsigned int iter_idx;

	vec_t* points; // Initial point + iterations
	vec_t* directions; // Initial direction + iterations

	unsigned int* presence;// [GRID_SIZE] ;
};
typedef struct data_t data_t;

struct intersection_t {
	vec_t pos, normal;
	double dist;
};
typedef struct intersection_t intersection_t;

/* FUNCTIONS DEF */

size_t coords2idx(const vec_t* coords);
int intersect_atom(vec_t* point, vec_t* dir, const vec_t* atom_pos, intersection_t* intersection);
void iter(data_t* data);

/* MAIN */

int main(int argc, char* argv[]) {
	size_t ITER_COUNT = 50;
	bool out_presence = false;
	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-iter") == 0) {
			if (argc < i + 1) {
				printf("-iter argument expected iteration count\n");
				return EXIT_FAILURE;
			}

			char* end;
			ITER_COUNT = strtoul(argv[i + 1], &end, 10);
			printf("ITER_COUNT = %llu\n", ITER_COUNT);

			++i;
		}
		else if (strcmp(argv[i], "-presence") == 0) {
			out_presence = true;
			printf("OUT_DATA = presence\n");
		}
		else {
			printf("Unrecognized option '%s'\n", argv[i]);
		}
	}
	
	data_t data = { .iter_idx = 0 };

	data.points = malloc((ITER_COUNT + 1) * sizeof(vec_t));
	if (!data.points) {
		printf("Failed to allocate %llu bytes\n", (ITER_COUNT + 1) * sizeof(vec_t));
		return EXIT_FAILURE;
	}

	data.directions = malloc((ITER_COUNT + 1) * sizeof(vec_t));
	if (!data.directions) {
		printf("Failed to allocate %llu bytes\n", (ITER_COUNT + 1) * sizeof(vec_t));
		free(data.points);
		return EXIT_FAILURE;
	}

	data.presence = malloc(GRID_SIZE * sizeof(unsigned int));
	if (!data.directions) {
		printf("Failed to allocate %llu bytes\n", GRID_SIZE * sizeof(unsigned int));
		free(data.points);
		free(data.directions);
		return EXIT_FAILURE;
	}
	for (size_t i = 0; i < GRID_SIZE; ++i)
		data.presence[i] = 0;
	
	/* Initial pos */
	data.points[0].x = 0;
	data.points[0].y = 0;
	/* Initial direction */
	data.directions[0].x = 0;
	data.directions[0].y = 1;
	vec_normalize(&data.directions[0]);

	for (size_t i = 0; i < ITER_COUNT; ++i)
		iter(&data);

	FILE* file = fopen("out.dat", "w");
	if (!file) {
		fprintf(stderr, "Failed to open 'out.dat'\n");

		free(data.points);
		free(data.directions);
		free(data.presence);
		return EXIT_FAILURE;
	}

	if (out_presence) {
		for (size_t i = 0; i < GRID_SIZE; ++i) {
			fprintf(file, "%u\n", data.presence[i]);
		}
	}
	else {
		for (size_t i = 0; i < ITER_COUNT + 1; ++i) {
			fprintf(file, "%f;%f\n", data.points[i].x, data.points[i].y);
		}
	}
	fclose(file);

	free(data.points);
	free(data.directions);
	free(data.presence);
	return EXIT_SUCCESS;
}

/* FUNCTIONS IMPL */

size_t coords2idx(const vec_t* coords) {
	const size_t x_idx = (size_t)floor((coords->x + 0.5) / GRID_INTERVAL);
	const size_t y_idx = (size_t)floor((coords->y + 0.5) / GRID_INTERVAL);

	const size_t grid_width  = (size_t)(1.0 / GRID_INTERVAL); // Grid of width 1 with intervals of 0.01

	return y_idx * grid_width + x_idx;
}

int intersect_atom(vec_t* point, vec_t* dir, const vec_t* atom_pos, intersection_t* intersection) {
	vec_t dist = vec_from_points(atom_pos, point);
	const double dot = vec_dot(dir, &dist);
	const double discr = 4.0 * (dot * dot - vec_length_sq(&dist) + (2.0 * ATOM_RADIUS) * (2.0 * ATOM_RADIUS)); // Twice the atom radius because there are fixed atom + moving one

	if (discr < 0)
		return 0;
	else {
		double dist = -0.5 * (2.0 * dot + sqrt(discr));
		if (dist <= 0) // If the target is behind the point's direction, we didnt actually hit
			return 0;

		intersection->dist = dist;

		const vec_t intersec = { .x = point->x + dir->x * dist, .y = point->y + dir->y * dist };
		intersection->pos.x = intersec.x;
		intersection->pos.y = intersec.y;

		intersection->normal = vec_from_points(atom_pos, &intersec);
		vec_normalize(&intersection->normal);
		return 1;
	}
}

void iter(data_t* data) {
	const vec_t atoms[4] = { {.x = -0.5, .y = 0.5 },
							 {.x = 0.5,  .y = 0.5 },
							 {.x = -0.5, .y = -0.5  },
							 {.x = 0.5,  .y = -0.5  } };

	intersection_t min_intersec = { .dist = -1.0 };

	for (size_t i = 0; i < 4; ++i) {
		intersection_t intersec;
		const int intersected = intersect_atom(&data->points[data->iter_idx], &data->directions[data->iter_idx], &atoms[i], &intersec);

		if (intersected && (min_intersec.dist < 0 || intersec.dist < min_intersec.dist))
			min_intersec = intersec;
	}

	vec_t tangent = { .x = min_intersec.normal.y, .y = -min_intersec.normal.x };
	double dir_dot_tan = vec_dot(&data->directions[data->iter_idx], &tangent);
	if (dir_dot_tan < 0) {
		tangent.x *= -1.0;
		tangent.y *= -1.0;
		dir_dot_tan *= -1.0;
	}

	vec_t presence_point = data->points[data->iter_idx];
	while (vec_dist_sq(&data->points[data->iter_idx], &presence_point) < min_intersec.dist * min_intersec.dist) {
		const size_t coords = coords2idx(&presence_point);
		data->presence[coords]++;

		presence_point.x += GRID_INTERVAL * data->directions[data->iter_idx].x;
		presence_point.y += GRID_INTERVAL * data->directions[data->iter_idx].y;
	}

	++data->iter_idx;

	data->points[data->iter_idx] = min_intersec.pos;
	data->directions[data->iter_idx].x = -(data->directions[data->iter_idx - 1].x - 2.0 * dir_dot_tan * tangent.x);
	data->directions[data->iter_idx].y = -(data->directions[data->iter_idx - 1].y - 2.0 * dir_dot_tan * tangent.y);
	vec_normalize(&data->directions[data->iter_idx]);
}