/****
 * @date Created on 2021-04-06 at 12:53:25 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing utilities to manipulate real double-precision vectors.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "common_util.h"

/**
 * Copies the content of the source vector, "src", to the destination vector "dst".
 * The two vectors are of the same length "n". The destination vector, "dst", is overwritten.
 */
void copy_vector(int n, double* src, double* dst) {
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		dst[i] = src[i];
	}
}

/**
 * Returns the Euclidean norm of the given array "x" of length "n" containing double precision values.
 */
double norm(int n, double* x) {
	double r2 = 0.;
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		r2 += x[i]*x[i];
	}
	return sqrt(r2);
}

/**
 * Computes the Euclidean distance between the two given vectors "x" and "y" of size "n".
 */
double distance(int n, double* x, double* y) {
	double dx, r2 = 0.;
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		dx = x[i] - y[i];
		r2 += dx*dx;
	}
	return sqrt(r2);
}

/**
 * Adds the vector "alpha*dx" to the vector "x". The last vector "x" is then modified.
 * The vectors are of length "n".
 */
void add_vector(int n, double* x, double* dx, double alpha) {
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		x[i] += alpha*dx[i];
	}
}

/**
 * Swaps the two components "i" and "j" of the given vector "x".
 */
void swap_vector(double* x, int i, int j) {
	double tmp = x[i];
	x[i] = x[j];
	x[j] = tmp;
}

/**
 * Sorts the given size-n real vector "x" in increasing order.
 */
void sort_vector(int n, double* x) {
	if (n <= 1) return; //Do not process too small arrays.
	int cmp(const void* a, const void* b) {//Comparison function.
		return *(double*)a > *(double*)b;
	};
	qsort(x, n, sizeof(double), cmp);
}

/**
 * Moves the duplicate entries of the size-n real vector "x" to the end of the array, and reduces its logical size "n".
 * This function does not change the relative order of the unique entries, and does not erase any entry of "x".
 * In addition, this function does not assume that the array "x" is sorted, and so uses a simple O(n^2) method.
 * @param "n" Logical size of the array "x".
 * @param "x" Array of double-precision real numbers.
 * @param "toler" Tolerance on the relative difference between any pair of entries of "x".
 * @returns "ndupl" The number of duplicate entries of "z" which have been removed.
 */
int uniq_vector(int* n, double* x, double toler) {
	if (*n <= 1) return 0; //Do not process too small arrays.
	int i, j, isuniq, size = 1;
	for (i = 1; i < *n; i++) {//Main loop.
		isuniq = 1;
		for (j = 0; j < size; j++) {//Loop on the sub-array to detect unicity.
			if (fabs(x[i] - x[j]) < toler*fabs(x[i]) || fabs(x[i] - x[j]) < toler) {//If the two numbers are close, then x[i] is not unique.
				isuniq = 0;
				break;
			}
		}
		if (isuniq) {//If "isuniq" is true, then save it into the sub-array.
			if (i != size) {//Swap only if necessary (in principle, size < i).
				//printf("[INFO] Swapping (%d,%d)...\n", i, size);
				swap_vector(x, i, size);
			}
			size++;
		}
	}
	int ndupl = *n - size;
	*n = size;
	return ndupl;
}

/**
 * Returns the q-quantile of the size-n real vector "x" using the linear interpolation method (without index excess).
 * This function assumes the input array "x" is already sorted in increasing order.
 */
double quantile(int n, double* x, double q) {
	if (n <= 0) {//Invalid array.
		return NAN;
	}
	else if (n == 1 || q <= 0.) {//If only one sample, then returns it.
		return x[0];
	}
	else if (q >= 1.) {//If q=1, then returns the last element.
		return x[n-1];
	}
	double idx = (n-1)*q;     //Continuous index in the interval [0, n-1].
	int i = (int)floor(idx);  //Integer part of the continuous index.
	return x[i] + (idx - i)*(x[i+1] - x[i]);
}

/**
 * Returns the mean of the given size-n real array "x".
 */
double mean_vector(int n, double* x) {
	if (n <= 0) return NAN; //Invalid array.
	double sum = 0.;
	int i;
	for (i = 0; i < n; i++) {
		sum += x[i];
	}
	return sum/n;
}

/**
 * Prints the given vector "x" of size "n" with name "name".
 */
void print_vector(int n, double* x, const char* name) {
	printf("[INFO] %s = {", name);
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		printf("%g", x[i]);
		if (i != n-1)
			printf(", ");
	}
	printf("}\n");
}

/**
 * Save the given vector "x" of size "n" with name "name" to the file "fp".
 */
void save_real_data(int n, double* x, const char* name, FILE* fp) {
	fprintf(fp, "%s=", name);
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		fprintf(fp, "%.16lg ", x[i]);
	}
	fprintf(fp, "\n");
}

/**
 * Counts the real numbers in the given string "arg".
 * The numbers are expected to end with a non-space character, followed by a space-like character or the null character.
 */
int count_real_data(const char* arg) {
	int i = 0;  //Character index.
	int n = 0;  //Current count of complex numbers.
	while (arg[i]) {//Loop on the string characters.
		if (!is_blank(arg[i]) && (is_blank(arg[i+1]) || arg[i+1] == '\0')) {
			n++;
		}
		i++;
	}
	return n;
}

/**
 * Parses the first "n" real numbers from the given data string "arg" to the array "data".
 * The expected format of "arg" is: "1.3 2.5 3. 1. -1.8 -4.7 +2.6", and so on.
 * Assumes that "data" is already allocated with enough place, i.e., at least "n" entries, and that the
 * string "arg" contains at least "n" valid numbers. This function also understands "nan" and "inf".
 */
void parse_real_data(int n, double* data, char* arg) {
	char *curarg = arg, *endarg;
	int i;
	for (i = 0; i < n; i++) {//Loop on the expected real numbers.
		data[i] = strtod(curarg, &endarg);
		if (endarg != NULL) {
			curarg = endarg;
		}
		else {//Progress to the next number.
			printf("[ERROR] Data corruption at '%.20s' (char %ld), aborting...\n", curarg, curarg-arg);
			exit(EXIT_FAILURE); //Continuing the program with corrupted data makes no sense.
		}
	}
}
