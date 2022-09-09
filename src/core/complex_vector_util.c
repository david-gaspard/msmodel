/****
 * @date Created on 2021-04-06 at 12:40:42 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to manipulate complex vectors.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "dcomplex_type.h"
#include "common_util.h"

/**
 * Copy the content of the given complex vector "src" of size "n" to the vector "dst".
 */
void copy_cvector(int n, dcomplex* src, dcomplex* dst) {
	int i;
	for (i = 0; i < n; i++) {
		dst[i] = src[i];
	}
}

/**
 * Returns the norm of the given complex vector "x" of size "n".
 */
double cnorm(int n, dcomplex* x) {
	double rex, imx, r2 = 0.;
	int i;
	for (i = 0; i < n; i++) {
		rex = creal(x[i]);
		imx = cimag(x[i]);
		r2 += rex*rex + imx*imx;
	}
	return sqrt(r2);
}

/**
 * Normalize the given complex vector "x" of size "n" to unit norm. Check whether the vector is zero in a safe way.
 * As a convention, the global complex phase of "x" is set such that the first nonzero element is a positive real number.
 * Returns the factor "fac" such that the vector "x" is equal to "fac" times the normalized "x".
 */
dcomplex cnormalize(int n, dcomplex* x) {
	dcomplex fac = 1.;
	int i = 0;
	while (i < n && x[i] == 0.) {//Find the first nonzero element.
		i++;
	}
	if (i < n) {//If i < n, then we found a nonzero element and the vector is not zero.
		fac = cnorm(n, x)*x[i]/cabs(x[i]);
		for (i = 0; i < n; i++) {
			x[i] /= fac;
		}
	}
	return fac;
}

/**
 * Computes the Euclidean distance between the two given complex vector "x" and "y" of size "n".
 */
double cdistance(int n, dcomplex* x, dcomplex* y) {
	double r2 = 0.;
	dcomplex dx;
	int i;
	for (i = 0; i < n; i++) {//Loop on the vector components.
		dx = x[i] - y[i];
		r2 += creal(dx)*creal(dx) + cimag(dx)*cimag(dx);
	}
	return sqrt(r2);
}

/**
 * Initialize the given vector "x" of size "n" with constant entries, and normalize it to unit norm.
 */
void constant_unit_cvector(int n, dcomplex* x) {
	int i;
	for (i = 0; i < n; i++) {
		x[i] = 1.;
	}
	cnormalize(n, x);
}

/**
 * Returns the dot product, "conj(x)·y", between the two given vectors "x" and "y" of size "n".
 * Note that the vector "x" is conjugated.
 */
dcomplex cdot_product(int n, dcomplex* x, dcomplex* y) {
	dcomplex res = 0.;
	int i;
	for (i = 0; i < n; i++) {
		res += conj(x[i])*y[i];
	}
	return res;
}

/**
 * Swaps the two complex numbers of given indices "i" and "j" in the array "z".
 * Useful for many array manipulation routines.
 */
void swap_cvector(dcomplex* z, int i, int j) {
	dcomplex tmp = z[i];
	z[i] = z[j];
	z[j] = tmp;
}

/**
 * Sorts the given size-n complex vector "z" in increasing order of the function "f(z)".
 */
void sort_cvector(int n, dcomplex* z, double (*f)(dcomplex)) {
	if (n <= 1) return; //Do not process too small arrays.
	int cmp(const void* a, const void* b) {//Comparison function.
		return f(*(dcomplex*)a) > f(*(dcomplex*)b);
	};
	qsort(z, n, sizeof(dcomplex), cmp);
}

/**
 * Moves the duplicate entries of the size-n complex vector "z" to the end of the array, and reduces its logical size "n".
 * This function does not change the relative order of the unique entries, and does not erase any entry of "z".
 * In addition, this function does not assume that the array "z" is sorted, and so uses a simple O(n^2) method.
 * @param "n" Logical size of the array "z".
 * @param "z" Array of double-precision complex numbers.
 * @param "toler" Tolerance on the relative difference between any pair of entries of "z".
 * @returns "ndupl" The number of duplicate entries of "z" which have been removed.
 */
int uniq_cvector(int* n, dcomplex* z, double toler) {
	if (*n <= 1) return 0; //Do not process too small arrays.
	int i, j, isuniq, size = 1;
	for (i = 1; i < *n; i++) {//Main loop.
		isuniq = 1;
		for (j = 0; j < size; j++) {//Loop on the sub-array to detect unicity.
			if (cabs(z[i] - z[j]) < toler*cabs(z[i]) || cabs(z[i] - z[j]) < toler) {//If the two numbers are close, then z[i] already exists in the sub-array.
				isuniq = 0;
				break;
			}
		}
		if (isuniq) {//If "isuniq" is true, then save it into the sub-array.
			if (i != size) {//Swap only if necessary (in principle, size < i).
				//printf("[INFO] Swapping (%d,%d)...\n", i, size);
				swap_cvector(z, i, size);
			}
			size++;
		}
	}
	int ndupl = *n - size;
	*n = size;
	return ndupl;
}

/**
 * Filters the entries of the given size-n complex vector "z" according to the flags of validity stored in the size-n "flag".
 * For every index "i", if flag[i] is equal to 1, then z[i] is saved. Otherwise, z[i] is moved to the end of the array "z" and its logical size "n" is reduced.
 * This function does not change the relative order of the valid entries. Note that the flag array is not re-ordered in consequence.
 * @param "n" Logical size of both arrays "z" and "flag". On exit, the reduced size of "z"
 * @param "z" Array of double-precision complex numbers to be filtered out. On exit, the filtered array.
 * @param "flag" Array of integers indicating whether or not the corresponding entry of "z" should be removed.
 * @returns "nfail" The number of failed entries in the array "z" which have been removed.
 */
int filter_cvector(int* n, dcomplex* z, int* flag) {
	int i, size = 0;
	for (i = 0; i < *n; i++) {//Loop on all the elements of the array.
		if (flag[i] == 1) {//If the i-th entry is valid.
			if (i != size) {//Only swap if necessary.
				//printf("[INFO] Swapping (%d,%d)...\n", i, size);
				swap_cvector(z, i, size); //Then add z[i] to the logical list, at the end.
			}
			size++;
		}
	}
	int nfail = *n - size;
	*n = size; //Saves the new size of the array of roots.
	return nfail;
}

/**
 * Prints the given vector "x" of size "n" with name "name".
 */
void print_cvector(int n, dcomplex* x, const char* name) {
	printf("[INFO] %s = {", name);
	int i;
	for (i = 0; i < n; i++) {
		printf("%g%+gi", creal(x[i]), cimag(x[i]));
		if (i != n-1)
			printf(", ");
	}
	printf("}\n");
}

/**
 * Returns the count of complex numbers in the given string "arg".
 * The numbers are expected to end with 'i' followed by a space-like character or the null character.
 */
int count_complex_data(const char* arg) {
	int i = 0;  //Character index.
	int n = 0;  //Current count of complex numbers.
	while (arg[i]) {
		if (arg[i] == 'i' && (is_blank(arg[i+1]) || arg[i+1] == '\0')) {
			n++;
		}
		i++;
	}
	return n;
}

/**
 * Parses the first "n" complex numbers from the given data string "arg" to the array "data".
 * The expected format of "arg" is: "1.3+2.5i 3+i -1.8-4.7i 2.6+5.4i 10.3+12.5i -5.1+4.2i", and so on.
 * Assumes that "data" is already allocated with enough place, i.e., at least "n" entries, and that the
 * string "arg" contains at least "n" valid numbers. This function also understands "nan" and "inf".
 */
void parse_complex_data(int n, dcomplex* data, char* arg) {
	char *curarg = arg, *endreal, *endimag;
	double re, im;
	int i;
	for (i = 0; i < n; i++) {//Loop on the expected complex numbers.
		re = strtod(curarg,  &endreal);
		im = strtod(endreal, &endimag);
		data[i] = re + im*I;
		if (endreal != NULL && is_char_in_str(endreal[0], "+- \t") && endimag != NULL && endimag[0] == 'i') {
			curarg = endimag + 1;
		}
		else {//Progress to the next number.
			printf("[ERROR] Data corruption at '%.20s' (char %ld), aborting...\n", curarg, curarg-arg);
			exit(EXIT_FAILURE); //Continuing the program with corrupted data makes no sense.
		}
	}
}
