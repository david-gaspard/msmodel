/****
 * @date Created on 2021-05-09 at 19:51:15 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing some assertion functions used for tests.
 ***/
#include <stdio.h>
#include <math.h>
#include "tag_strings.h"
#include "real_vector_util.h"
#include "complex_vector_util.h"

/**
 * Prints to the console the name of the test, or the name of the function which is tested.
 */
void print_test(const char* name) {
	printf("\033[1m====== TEST: %s ======\033[0m\n", name);
}

/**
 * Checks if the two real vectors are equal. Returns 1 on failure, 0 otherwise.
 */
int assert_vector_equals(int n, double* x, int n_expc, double* x_expc) {
	int i, failed = 0;
	if (n != n_expc) {
		failed = 1;
		printf(STR_FAIL "Invalid size=%d, expected=%d \n", n, n_expc);
		if (n_expc < n) //take the minimum of both to avoid out of bound error.
			n = n_expc;
	}
	for (i = 0; i < n; i++) {
		if (x[i] != x_expc[i]) {
			failed = 1;
			printf(STR_FAIL "Invalid element x[%d]=%g, expected=%g.\n", i, x[i], x_expc[i]);
			print_vector(n_expc, x_expc, "Expected x");
			break;
		}
	}
	if (failed == 0) {
		printf(STR_PASS "Vector is correct.\n");
	}
	return failed;
}

/**
 * Checks if the two complex vectors are equal. Returns 1 on failure, 0 otherwise.
 */
int assert_cvector_equals(int n, dcomplex* z, int n_expc, dcomplex* z_expc) {
	int i, failed = 0;
	if (n != n_expc) {
		failed = 1;
		printf(STR_FAIL "Invalid size=%d, expected=%d \n", n, n_expc);
		if (n_expc < n) //take the minimum of both to avoid out of bound error.
			n = n_expc;
	}
	for (i = 0; i < n; i++) {
		if (z[i] != z_expc[i]) {
			failed = 1;
			printf(STR_FAIL "Invalid element z[%d]=%g%+gi, expected=%g%+gi.\n", i, creal(z[i]), cimag(z[i]), creal(z_expc[i]), cimag(z_expc[i]));
			print_cvector(n_expc, z_expc, "Expected z");
			break;
		}
	}
	if (failed == 0) {
		printf(STR_PASS "Vector is correct.\n");
	}
	return failed;
}

/**
 * Check that the real vector "x" is almost equal to "x_expc" to within the global tolerance.
 */
int assert_vector_close(int n, double* x, int n_expc, double* x_expc, double toler) {
	int i, failed = 0;
	if (n != n_expc) {
		failed = 1;
		printf(STR_FAIL "Invalid size=%d, expected=%d \n", n, n_expc);
		if (n_expc < n) //take the minimum of both to avoid out of bound error.
			n = n_expc;
	}
	double err, abstol;
	for (i = 0; i < n; i++) {
		err = fabs(x[i] - x_expc[i]);
		abstol = toler*fabs(x[i]); //Absolute tolerance.
		if (err > abstol) {
			failed = 1;
			printf(STR_FAIL "Invalid element x[%d]=%.15g, expected=%.15g, err=%g (abstol=%g).\n", i, x[i], x_expc[i], err, abstol);
			print_vector(n_expc, x_expc, "Expected x");
			break;
		}
	}
	if (failed == 0) {
		printf(STR_PASS "Vector is correct.\n");
	}
	return failed;
}

/**
 * Check that the complex vector "z" is almost equal to "z_expc" to within the global tolerance.
 */
int assert_cvector_close(int n, dcomplex* z, int n_expc, dcomplex* z_expc, double toler) {
	int i, failed = 0;
	if (n != n_expc) {
		failed = 1;
		printf(STR_FAIL "Invalid size=%d, expected=%d \n", n, n_expc);
		if (n_expc < n) //take the minimum of both to avoid out of bound error.
			n = n_expc;
	}
	double err, abstol;
	for (i = 0; i < n; i++) {
		err = cabs(z[i] - z_expc[i]);
		abstol = toler*cabs(z[i]); //Absolute tolerance.
		if (err > abstol) {
			failed = 1;
			printf(STR_FAIL "Invalid element z[%d]=%.15g%+.15gi, expected=%.15g%+.15gi, err=%g (abstol=%g).\n",
				i, creal(z[i]), cimag(z[i]), creal(z_expc[i]), cimag(z_expc[i]), err, abstol);
			print_cvector(n_expc, z_expc, "Expected z");
			break;
		}
	}
	if (failed == 0) {
		printf(STR_PASS "Vector is correct.\n");
	}
	return failed;
}

/**
 * Checks if the size-n array "z" is a subset of "z_expc". Returns 1 on failure, 0 otherwise.
 * There is a failure when at least one element of "z" does not belong to "z_expc".
 */
int assert_cvector_subset(int n, dcomplex* z, int n_expc, dcomplex* z_expc, double toler) {
	int i, j, jmatch, failed = 0;
	double minerr, err, abstol;
	for (i = 0; i < n; i++) {//Loop on the elements of "z".
		jmatch = 0; //Matching element in "z_expc".
		minerr = cabs(z[i] - z_expc[0]);
		abstol = toler*cabs(z[i]); //Absolute tolerance.
		for (j = 1; j < n_expc; j++) {//Loop on the expected elements to find the corresponding one.
			err = cabs(z[i] - z_expc[j]);
			if (err < minerr) {//Find the minimum error.
				minerr = err;
				jmatch = j;
			}
		}
		if (minerr > abstol) {//If there is no satisfactory match, then failed.
			failed = 1;
			printf(STR_FAIL "Invalid element z[%d]=%.15g%+.15gi, expected=%.15g%+.15gi, err=%g (abstol=%g).\n",
				i, creal(z[i]), cimag(z[i]), creal(z_expc[jmatch]), cimag(z_expc[jmatch]), minerr, abstol);
			print_cvector(n_expc, z_expc, "Expected");
			break;
		}
	}
	if (failed == 0) {
		printf(STR_PASS "All elements belong to the expected array.\n");
	}
	return failed;
}
