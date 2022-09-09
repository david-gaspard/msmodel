/****
 * @date Created on 2021-03-12 at 16:49:14 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the complex domain utilities.
 ***/
//#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "domain_type.h"
#include "common_util.h"

/**
 * Initializes the domain "dom" from the given source domain "src".
 * Also check the validity of the domain boundaries. Returns 1 on failure, i.e., invalid domain boundaries.
 */
int init_domain(Domain* dom, double xmin, double xmax, double ymin, double ymax) {
	int fail = 0;
	if (xmin >= xmax) {
		printf("[ERROR] Invalid domain x-range, xmin=%g and xmax=%g...\n", xmin, xmax);
		fail = 1;
	}
	if (ymin >= ymax) {
		printf("[ERROR] Invalid domain y-range, ymin=%g and ymax=%g...\n", ymin, ymax);
		fail = 1;
	}
	dom->xmin = xmin;
	dom->xmax = xmax;
	dom->ymin = ymin;
	dom->ymax = ymax;
	return fail;
}

/**
 * Copy the source domain "src" to the destination "dst".
 */
void copy_domain(Domain* src, Domain* dst) {
	init_domain(dst, src->xmin, src->xmax, src->ymin, src->ymax);
}

/**
 * Parses the domain "dom" from the "narg" arguments contained in "args".
 * The "dom" pointer is supposed to be already allocated and can be dereferenced.
 * Returns 0 on success (successfully parsed xrange/yrange) and 1 on failure (no xrange/yrange found).
 */
int parse_domain(Domain* dom, int narg, char** args) {
	int fail = 0; //The xrange/yrange pair is found by default.
	double xmin = 0., xmax = 0., ymin = 0., ymax = 0.;
	char* value = get_value(narg, args, "xrange");
	if (value[0]) {
		if (sscanf(value, "%lg:%lg", &xmin, &xmax) != 2) {
			printf("[ERROR] Invalid domain x-range '%s'...\n", value);
			fail = 1;
		}
	}
	else {//The xrange is not found.
		//printf("[ERROR] Domain x-range not found...\n");
		fail = 1;
	}
	value = get_value(narg, args, "yrange");
	if (value[0]) {
		if (sscanf(value, "%lg:%lg", &ymin, &ymax) != 2) {
			printf("[ERROR] Invalid domain y-range '%s'...\n", value);
			fail = 1;
		}
	}
	else {//The yrange is not found.
		//printf("[ERROR] Domain y-range not found...\n");
		fail = 1;
	}
	fail |= init_domain(dom, xmin, xmax, ymin, ymax);
	return fail;
}

/**
 * Saves the given domain "dom" to the given file.
 */
void save_domain(Domain* dom, FILE* fp) {
	fprintf(fp, "xrange=%.16lg:%.16lg\n", dom->xmin, dom->xmax);
	fprintf(fp, "yrange=%.16lg:%.16lg\n", dom->ymin, dom->ymax);
}

/**
 * Returns the length of the domain "dom" along the real direction.
 */
double length_x(Domain* dom) {
	return fabs(dom->xmax - dom->xmin);
}

/**
 * Returns the length of the domain "dom" along the imaginary direction.
 */
double length_y(Domain* dom) {
	return fabs(dom->ymax - dom->ymin);
}

/**
 * Returns the surface area of the given domain "dom".
 */
double domain_area(Domain* dom) {
	return length_x(dom)*length_y(dom);
}

/**
 * Returns 1 if the given 2D point "p" lies in the complex domain "dom".
 */
int is_point_in_domain(Domain* dom, double* p) {
	return dom->xmin < p[0] && p[0] < dom->xmax && dom->ymin < p[1] && p[1] < dom->ymax;
}

/**
 * Returns 1 if the given complex number "z" lies in the complex domain "dom".
 */
int is_in_domain(Domain* dom, dcomplex z) {
	return dom->xmin < creal(z) && creal(z) < dom->xmax && dom->ymin < cimag(z) && cimag(z) < dom->ymax;
}

/**
 * Moves the complex numbers from the array "z", which are located outside of the domain "dom", to the end of the array "z".
 * The array length "n" is then reduced accordingly. This leads to a smaller array but without removing anything,
 * in other words the logical size of the array "z" becomes smaller than its actual capacity.
 * This routine does not change the relative order of the valid elements in "z".
 * @param "dom" Domain in which the complex points of "z" are expected to be.
 * @param "n" Size of the complex array "z" at the beginning. On exit, number of values of "z" in the domain "dom".
 * @param "z" On exit, only points inside the domain are kept in the size-n array.
 * @param "nexcl" The number of points removed from the array "z".
 */
int exclude_from_domain(Domain* dom, int* n, dcomplex* z) {
	int i, size = 0;  //Logical size.
	dcomplex tmp;
	for (i = 0; i < *n; i++) {//Loop on all the elements of the array.
		if (is_in_domain(dom, z[i])) {//If z[i] is inside the domain.
			if (i != size) {//Swap only if necessary.
				tmp = z[i]; //Then add z[i] to the logical list, at the end.
				z[i] = z[size];
				z[size] = tmp;
			}
			size++;
		}
	}
	int nexcl = *n - size;
	*n = size;  //Reduce the logical array length.
	return nexcl;
}

/**
 * Gives the i-th Halton-van der Corput number in base "b" equi-distributed in the unit interval [0, 1].
 * This can be used to generate the (2,3) Halton points sequence for instance.
 * Based on: https://en.wikipedia.org/wiki/Van_der_Corput_sequence, and similar articles.
 */
static double halton_number(int i, int b) {//Private function
	double h = 0., f = 1.;
	while (i > 0) {
		f /= b;
		h += f*(i%b);
		i /= b;
	}
	return h;
}

/**
 * Returns the i-th (2,3) Halton point in the given complex domain "dom".
 * Note that this sequence usually begins at i=1, not at i=0. Indeed, i=0 is the corner point (xmin, ymin).
 */
dcomplex halton_point_2_3(Domain* dom, int i) {
	return dom->xmin + (dom->xmax - dom->xmin)*halton_number(i, 2) + I*(dom->ymin + (dom->ymax - dom->ymin)*halton_number(i, 3));
}

/**
 * Stores the "nz" first Halton points of the given domain "dom" into the given array "z".
 * @param "dom" Domain on which the Halton points are computed.
 * @param "n" Number of desired Halton points.
 * @param "z" On exit, array containing the halton points (size "n").
 */
void set_halton_points(Domain* dom, int n, dcomplex* z) {
	int i;
	for (i = 0; i < n; i++) {
		z[i] = halton_point_2_3(dom, i+1); //Pick up sample points inthe domain "dom" using low-discrepancy Halton-van der Corput's (2,3) sequence.
	}
}
