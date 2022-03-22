#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

#include "vec.h"

#define ATOM_RADIUS 0.25 // As a proportion of lattice's size
#define GRID_SIZE 200

/* STRUCTS */
struct data_t {	
	unsigned int iter_idx;

	vec_t* points; // Initial point + iterations
	vec_t* directions; // Initial direction + iterations

	unsigned int* presence; // [GRID_SIZE^3];
	
	double* thetas;
	double* dots;

	const vec_t atoms[4];
};
typedef struct data_t data_t;

struct intersection_t {
	vec_t pos, normal;
	double dist;
};
typedef struct intersection_t intersection_t;

/* FUNCTIONS DEF */

size_t coords2idx(const vec_t* coords);
int init_data(data_t* data, size_t ITER_COUNT);

int intersect_atom(vec_t* point, vec_t* dir, const vec_t* atom_pos, intersection_t* intersection);
int iter_wall(data_t* data);

double potential(vec_t* pos, const vec_t* atoms_pos);
double dpot_dx(vec_t* pos, const vec_t* atoms_pos);
double dpot_dy(vec_t* pos, const vec_t* atoms_pos);
int iter_pot(data_t* data);

/* MAIN */

int main(int argc, char* argv[]) {
	size_t ITER_COUNT = 50;
	bool TRUE_POTENTIAL = false;
	double x_i = 0., y_i = 0., dx_i = 1., dy_i = 1.;

	for (int i = 1; i < argc; ++i) {
		if (strcmp(argv[i], "-iter") == 0) {
			if (argc < i + 1) {
				printf("-iter argument expected iteration count\n");
				return EXIT_FAILURE;
			}

			ITER_COUNT = strtoul(argv[i + 1], NULL, 10);

			++i;
		}
		else if (strcmp(argv[i], "-init") == 0) {
			if (argc < i + 4) {
				printf("-init argument expected initial position and direction\n");
				return EXIT_FAILURE;
			}

			x_i  = strtod(argv[i + 1], NULL);
			y_i  = strtod(argv[i + 2], NULL);
			dx_i = strtod(argv[i + 3], NULL);
			dy_i = strtod(argv[i + 4], NULL);

			i += 4;
		}
		else if (strcmp(argv[i], "-pot") == 0) {
			TRUE_POTENTIAL = true;
		}
		else {
			printf("Unrecognized option '%s'\n", argv[i]);
		}
	}

	printf("ITER_COUNT = %llu\n", ITER_COUNT);
	printf("TRUE_POTENTIAL = %d\n", TRUE_POTENTIAL);

	data_t data = {
		.iter_idx = 0,
		.atoms = { {.x = -0.5, .y = 0.5   },
		           {.x = 0.5,  .y = 0.5   },
		           {.x = -0.5, .y = -0.5  },
		           {.x = 0.5,  .y = -0.5  } },
	};

	if (init_data(&data, ITER_COUNT) == -1) {
		return EXIT_FAILURE;
	}

	/* Initial pos */
	data.points[0].x = x_i;
	data.points[0].y = y_i;
	/* Initial direction */
	data.directions[0].x = dx_i;
	data.directions[0].y = dy_i;
	vec_normalize(&data.directions[0]);
	
	printf("(x_i, y_i) = (%f, %f)\n(dx_i, dy_i) = (%f, %f)\n", x_i, y_i, data.directions[0].x, data.directions[0].y);

	for (size_t i = 0; i < ITER_COUNT; ++i) {
		int iter = 0;
		if (TRUE_POTENTIAL)
			iter = iter_pot(&data);
		else
			iter = iter_wall(&data);

		if (iter == -1) {
			free(data.points);
			free(data.directions);
			free(data.presence);
			free(data.thetas);
			free(data.dots);
			return EXIT_FAILURE;
		}
	}

	FILE* file = fopen("presence.dat", "w");
	if (!file) {
		fprintf(stderr, "Failed to open 'presence.dat'\n");

		free(data.points);
		free(data.directions);
		free(data.presence);
		free(data.thetas);
		free(data.dots);
		return EXIT_FAILURE;
	}
	for (size_t i = 0; i < GRID_SIZE * GRID_SIZE; ++i) {
		fprintf(file, "%u\n", data.presence[i]);
	}
	fclose(file);

	file = fopen("trajectories.dat", "w");
	if (!file) {
		fprintf(stderr, "Failed to open 'trajectories.dat'\n");

		free(data.points);
		free(data.directions);
		free(data.presence);
		free(data.thetas);
		free(data.dots);
		return EXIT_FAILURE;
	}
	for (size_t i = 0; i < ITER_COUNT + 1; ++i) {
		fprintf(file, "%f;%f\n", data.points[i].x, data.points[i].y);
	}
	fclose(file);

	file = fopen("phase.dat", "w");
	if (!file) {
		fprintf(stderr, "Failed to open 'trajectories.dat'\n");

		free(data.points);
		free(data.directions);
		free(data.presence);
		free(data.thetas);
		free(data.dots);
		return EXIT_FAILURE;
	}
	for (size_t i = 0; i < ITER_COUNT; ++i) {
		fprintf(file, "%f;%f\n", data.thetas[i], data.dots[i]);
	}
	fclose(file);

	free(data.points);
	free(data.directions);
	free(data.presence);
	free(data.thetas);
	free(data.dots);
	return EXIT_SUCCESS;
}

/* FUNCTIONS IMPL */

size_t coords2idx(const vec_t* coords) {
	const size_t x_idx = (size_t)floor(fmax(0.0, fmin(1.0, coords->x + 0.5)) * (GRID_SIZE - 1.0));
	const size_t y_idx = (size_t)floor(fmax(0.0, fmin(1.0, coords->y + 0.5)) * (GRID_SIZE - 1.0));

	return y_idx * GRID_SIZE + x_idx;
}

int init_data(data_t* data, size_t ITER_COUNT) {
	data->points = malloc((ITER_COUNT + 1) * sizeof(vec_t));
	if (!data->points) {
		printf("Failed to allocate %llu bytes\n", (ITER_COUNT + 1) * sizeof(vec_t));
		return -1;
	}

	data->directions = malloc((ITER_COUNT + 1) * sizeof(vec_t));
	if (!data->directions) {
		printf("Failed to allocate %llu bytes\n", (ITER_COUNT + 1) * sizeof(vec_t));
		free(data->points);
		return -1;
	}

	data->presence = malloc(GRID_SIZE * GRID_SIZE * sizeof(unsigned int));
	if (!data->directions) {
		printf("Failed to allocate %llu bytes\n", GRID_SIZE * GRID_SIZE * sizeof(unsigned int));
		free(data->points);
		free(data->directions);
		return -1;
	}
	for (size_t i = 0; i < GRID_SIZE * GRID_SIZE; ++i)
		data->presence[i] = 0;

	data->thetas = malloc(ITER_COUNT * sizeof(double));
	if (!data->thetas) {
		printf("Failed to allocate %llu bytes\n", ITER_COUNT * sizeof(double));

		free(data->points);
		free(data->directions);
		free(data->presence);
		return -1;
	}

	data->dots = malloc(ITER_COUNT * sizeof(double));
	if (!data->dots) {
		printf("Failed to allocate %llu bytes\n", ITER_COUNT * sizeof(double));

		free(data->points);
		free(data->directions);
		free(data->presence);
		free(data->dots);
		return -1;
	}

	return 0;
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

int iter_wall(data_t* data) {
	intersection_t min_intersec = { .dist = -1.0 };

	for (size_t i = 0; i < 4; ++i) {
		intersection_t intersec;
		const int intersected = intersect_atom(&data->points[data->iter_idx], &data->directions[data->iter_idx], &data->atoms[i], &intersec);

		if (intersected && (min_intersec.dist < 0 || intersec.dist < min_intersec.dist))
			min_intersec = intersec;
	}

	if (min_intersec.dist < 0) {
		printf("[ERROR] Point (%f,%f) in direction (%f,%f) never intersected!\n", data->points[data->iter_idx].x, data->points[data->iter_idx].y,
			                                                                      data->directions[data->iter_idx].x, data->directions[data->iter_idx].y);
		return -1;
	}

	vec_t presence_point = data->points[data->iter_idx];
	while (vec_dist_sq(&data->points[data->iter_idx], &presence_point) < min_intersec.dist * min_intersec.dist) {
		const size_t coords = coords2idx(&presence_point);
		data->presence[coords]++;

		presence_point.x += data->directions[data->iter_idx].x / GRID_SIZE;
		presence_point.y += data->directions[data->iter_idx].y / GRID_SIZE;
	}

	const vec_t tangent = { .x = min_intersec.normal.y, .y = -min_intersec.normal.x };
	const double dir_dot_tan = vec_dot(&data->directions[data->iter_idx], &tangent);

	++data->iter_idx;

	data->points[data->iter_idx] = min_intersec.pos;

	// Formula used for reflections: u' = -(u - 2(u . t)t)
	data->directions[data->iter_idx].x = -(data->directions[data->iter_idx - 1].x - 2.0 * dir_dot_tan * tangent.x);
	data->directions[data->iter_idx].y = -(data->directions[data->iter_idx - 1].y - 2.0 * dir_dot_tan * tangent.y);
	vec_normalize(&data->directions[data->iter_idx]);

	// θ = arccos(OP . (1,0))
	const vec_t unitary_x = { .x = 0, .y = 1 };
	vec_normalize(&min_intersec.pos);
	const vec_t OP = vec_from_points(&VEC_ZERO, &min_intersec.pos);
	double theta = acos(vec_dot(&OP, &unitary_x));

	if (min_intersec.pos.x > 0) {
		theta *= -1.;
	}

	data->thetas[data->iter_idx - 1] = theta;

	// (u . t)
	const double dot = vec_dot(&data->directions[data->iter_idx - 1], &tangent);
	data->dots[data->iter_idx - 1] = dot;

	return 0;
}

double potential(vec_t* pos, const vec_t* atoms_pos) {
	const double D = 0.4;
	const double a = 1.0 / (2.0 * ATOM_RADIUS);
	const double ATOMS_DIST = 1;

	double pot = 0.0;

	for (size_t i = 0; i < 4; ++i) {
		const double DIST = vec_dist(&atoms_pos[i], pos);
		pot += (1 - exp(-a * (DIST - ATOMS_DIST))) * (1 - exp(-a * (DIST - ATOMS_DIST)));
	}

	return D * pot;
}

double dpot_dx(vec_t* pos, const vec_t* atom_pos) {
	const double dx = 0.5 / GRID_SIZE;
	const double pot_x = potential(pos, atom_pos);
	pos->x += dx;
	
	const double dv = (potential(pos, atom_pos) - pot_x) / dx;
	pos->x -= dx;

	return dv;
}

double dpot_dy(vec_t* pos, const vec_t* atom_pos) {
	const double dy = 0.5 / GRID_SIZE;
	const double pot_y = potential(pos, atom_pos);
	pos->y += dy;

	const double dv = (potential(pos, atom_pos) - pot_y) / dy;
	pos->y -= dy;

	return dv;
}

int iter_pot(data_t* data) {
	const size_t i = data->iter_idx;

	const size_t coords = coords2idx(&data->points[i]);
	data->presence[coords]++;

	// Actually -strength
	const double strength_x = dpot_dx(&data->points[i], data->atoms);
	const double strength_y = dpot_dy(&data->points[i], data->atoms);

	// To roughly move from d²x = 1/GRID_SIZE each iteration when in null-potential space
	const double speed_mag = vec_length(&data->directions[i]);
	const double strength_mag = sqrt(strength_x * strength_x + strength_y * strength_y);

	const double dt = speed_mag / strength_mag * (sqrt(1 + 2 * strength_mag / (GRID_SIZE * speed_mag * speed_mag)) - 1);

	data->points[i + 1].x = data->points[i].x + data->directions[i].x * dt - 0.5 * strength_x * dt * dt;
	data->points[i + 1].y = data->points[i].y + data->directions[i].y * dt - 0.5 * strength_y * dt * dt;

	data->directions[i + 1].x = data->directions[i].x - 0.5 * (strength_x + dpot_dx(&data->points[i + 1], data->atoms)) * dt;
	data->directions[i + 1].y = data->directions[i].y - 0.5 * (strength_y + dpot_dy(&data->points[i + 1], data->atoms)) * dt;
	
	++data->iter_idx;

	return 0;
}