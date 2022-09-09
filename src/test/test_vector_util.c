/****
 * @date Created on 2021-03-28 at 12:27:27 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing in-depth tests for the complex and real vector tools, and especially data parsing functions.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tag_strings.h"
#include "assertion.h"
#include "real_vector_util.h"
#include "complex_vector_util.h"

#define TOLER      1.e-14     //General tolerance used for tests on relative errors of numerical computations.

/**************************************
 * TESTING THE REAL VECTORS UTILITIES:
 *************************************/
/**
 * Test the norm() function against some expected results.
 */
int test_norm_from_input(int n, double* x, double normx_expc) {
	print_vector(n, x, "Input x");
	double normx = norm(n, x);
	int failed = (normx != normx_expc);
	if (failed) {
		printf(STR_FAIL "Invalid norm=%g, expected=%g \n", normx, normx_expc);
	}
	else {
		printf(STR_PASS "Norm(x) = %g \n", normx);
	}
	return failed;
}
/**
 * Make all the tests of the norm() function.
 */
int test_norm() {
	print_test("norm()");
	double x1[] = {1, 2, 3, 4};
	double normx1_expc = sqrt(30.);
	double x2[] = {6.6, 4., 2.4, 3., 8.3, 4.4, 7.5, 9.9, -8.3, -0.8};
	double normx2_expc = sqrt(386.36);
	double x3[] = {0, 1, 0, 0};
	double normx3_expc = 1;
	
	int nf = 0; //Number of failed test.
	nf += test_norm_from_input(sizeof(x1)/sizeof(x1[0]), x1, normx1_expc);
	nf += test_norm_from_input(sizeof(x2)/sizeof(x2[0]), x2, normx2_expc);
	nf += test_norm_from_input(sizeof(x3)/sizeof(x3[0]), x3, normx3_expc);
	if (nf == 0) {
		printf(STR_PASS "norm(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the distance between two real vectors.
 */
int test_distance() {
	print_test("distance()");
	int nf = 0; //Number of failed test.
	double vec[] = {2, 0, -3, 7, -2, 4};
	int n = sizeof(vec)/(2*sizeof(vec[0])); //n=3.
	print_vector(n, vec,   "x");
	print_vector(n, vec+n, "y");
	double dist = distance(n, vec, vec+n);
	double dist_expc = sqrt(78.);
	printf("[INFO] Distance |x-y| = %g\n", dist);
	if (fabs(dist/dist_expc - 1.) > TOLER) {
		nf++;
		printf(STR_FAIL "Invalid distance=%g, expected=%g.\n", dist, dist_expc);
	}
	dist = distance(n, vec, vec);
	dist_expc = 0.;
	printf("[INFO] Distance |x-x| = %g\n", dist);
	if (dist > TOLER) {
		nf++;
		printf(STR_FAIL "Invalid distance=%g, expected=%g.\n", dist, dist_expc);
	}
	if (nf == 0) {
		printf(STR_PASS "distance(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the sorting of real vectors.
 */
int test_sort_vector() {
	print_test("sort_vector()");
	double x[] = {-5.1, 1.6, -4.5, 4.3, 7.9, -5.4, 6.3, 0.5, 8.2, -7.1, 11, -8.4, 2.3, -1.8};
	int n = sizeof(x)/sizeof(x[0]);
	print_vector(n, x, "Input  x");
	sort_vector(n, x);
	print_vector(n, x, "Sorted x");
	
	int i, failed = 0;
	for (i = 1; i < n; i++) {
		if (x[i-1] > x[i]) {//The array must be sorted in increasing order by convention.
			failed = 1;
			printf(STR_FAIL "Invalid order between x[%d]=%g and x[%d]=%g, expected x[i-1]<x[i].\n", i-1, x[i-1], i, x[i]);
		}
	}
	if (failed == 0) {
		printf(STR_PASS "sort_vector(): All tests passed.\n");
	}
	return failed;
}

/**
 * Test the deletion of duplicate entries in a real vector.
 */
int test_uniq_vector_from_input(int n, double* x, double toler, int n_expc, double* x_expc) {
	//int n_bak = n; //Backup the original number of components.
	print_vector(n, x, "Input x");
	uniq_vector(&n, x, toler);
	print_vector(n, x, "Uniqs x");
	//print_vector(n_bak, x, "Whole x");
	return assert_vector_equals(n, x, n_expc, x_expc);
}

/**
 * Make all the test on the uniq_vector() function.
 */
int test_uniq_vector() {
	print_test("uniq_vector()");
	double x1[] = {8, 8, 8, 1, 8, 10, 9, 9, 10, 10, 1, 9, 10, 10, 10, 9, 9, 1, 8, 8, 1, 1, 8};
	double x1_expc[] = {8, 1, 10, 9};
	double toler1 = 1e-15; //Tolerate nearly zero relative difference.
	
	double x2[] = {-1, 3.140, 3, 2, 3, 7, 7, -2.5, -1, 3, 2, 3.145, 7.001, 3};
	double x2_expc[] = {-1, 3.14, 3, 2, 7, -2.5};
	double toler2 = 0.01; //Tolerate 1% relative difference.
	
	double x3[] = {0, 0, 3, 3, 0, 3, 2, 2, 2, 3.01, 3, 2, 3, 2, 2.1, 2, 2, 3};
	double x3_expc[] = {0, 3, 2};
	double toler3 = 0.1; //Tolerate 10% relative difference.
	
	int nf = 0;
	nf += test_uniq_vector_from_input(sizeof(x1)/sizeof(x1[0]), x1, toler1, sizeof(x1_expc)/sizeof(x1_expc[0]), x1_expc);
	nf += test_uniq_vector_from_input(sizeof(x2)/sizeof(x2[0]), x2, toler2, sizeof(x2_expc)/sizeof(x2_expc[0]), x2_expc);
	nf += test_uniq_vector_from_input(sizeof(x3)/sizeof(x3[0]), x3, toler3, sizeof(x3_expc)/sizeof(x3_expc[0]), x3_expc);
	if (nf == 0) {
		printf(STR_PASS "uniq_vector(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the computation of the quantile in a single real array.
 */
int test_quantile_from_input(int n, double* x, double q, double expc_res) {
	//print_vector(n, x, "Input x");
	double res = quantile(n, x, q);
	if (fabs(1. - res/expc_res) > TOLER) {
		printf(STR_FAIL "Wrong %g-quantile %.16g, expected %.16g\n", q, res, expc_res);
		return 1;
	}
	else {
		printf(STR_PASS "Correct %g-quantile %.16g\n", q, res);
		return 0;
	}
}

/**
 * Test the quantile function on already sorted arrays.
 */
int test_quantile() {
	print_test("quantile()");
	double x1[] = {1, 2, 2, 3, 4, 5, 6, 7, 7, 7, 8, 9, 12, 13, 15, 20};
	int n1 = sizeof(x1)/sizeof(x1[0]);
	double x1_quartile1 = 3.75;
	double x1_quartile2 = 7.00; //median
	double x1_quartile3 = 9.75;
	
	double x2[] = {0.015045, 0.049862, 0.069981, 0.076965, 0.077408, 0.080334, 0.081385, 0.112491, 0.133678,
		0.15692, 0.157094, 0.17844, 0.214606, 0.218942, 0.227589, 0.239692, 0.245154, 0.275527, 0.300003,
		0.341011, 0.431336, 0.469651, 0.481844, 0.495878, 0.52273, 0.526701, 0.527372, 0.532473, 0.586434,
		0.64529, 0.646391, 0.701921, 0.708026, 0.719218, 0.744694, 0.774961, 0.777044, 0.777319, 0.810449,
		0.811952, 0.820832, 0.824305, 0.827071, 0.8327, 0.934049, 0.936541, 0.939792, 0.953401, 0.97247, 0.997978};
	int n2 = sizeof(x2)/sizeof(x2[0]);
	double x2_quartile1 = 0.21569;
	double x2_quartile2 = 0.5247155; //median
	double x2_quartile3 = 0.77725025;
	
	int nf = 0;
	nf += test_quantile_from_input(n1, x1, 0.25, x1_quartile1);
	nf += test_quantile_from_input(n1, x1, 0.50, x1_quartile2);
	nf += test_quantile_from_input(n1, x1, 0.75, x1_quartile3);
	nf += test_quantile_from_input(n2, x2, 0.25, x2_quartile1);
	nf += test_quantile_from_input(n2, x2, 0.50, x2_quartile2);
	nf += test_quantile_from_input(n2, x2, 0.75, x2_quartile3);
	if (nf == 0) {
		printf(STR_PASS "quantile(): All tests passed.\n");
	}
	return nf;
}

/**
 * Tests the output of "strtod" with the given "number".
 */
int test_strtod_from_input(double number) {
	int failed = 1; //This test is failed by default.
	char str[50];
	sprintf(str, "\t\n %g_this is test...", number);
	char *ptr;
	double result = strtod(str, &ptr);
	if (result != number && !isnan(number)) {
		printf(STR_FAIL "test_strtod: Result=%g is not equal to expected=%g...\n", result, number);
	}
	else if (ptr == NULL || ptr[0] != '_') {
		printf(STR_FAIL "test_strtod: String=|%s| is not expected...\n", ptr);
	}
	else {
		failed = 0;
		printf(STR_PASS "test_strtod(): Passed with %g.\n", number);
	}
	return failed;
}

/**
 * Test many possible inputs for strtod().
 */
int test_strtod() {
	print_test("strtod()");
	int nf = 0; //Number of failed tests.
	nf += test_strtod_from_input(+0.0);
	nf += test_strtod_from_input(-0.0);
	nf += test_strtod_from_input(10);
	nf += test_strtod_from_input(-0.00012458);
	nf += test_strtod_from_input(1.18e-7);
	nf += test_strtod_from_input(-1.245e10);
	nf += test_strtod_from_input(NAN);
	nf += test_strtod_from_input(INFINITY);
	nf += test_strtod_from_input(-INFINITY);
	return nf;
}

/**
 * Test the real array parsing function parse_real_data().
 */
int test_parse_real_data_from_input(char* arg, int n_expc, double* x_expc) {
	printf("[INFO] Input string = '%s'\n", arg);
	int n = count_real_data(arg); //Test the estimate of arguments number.
	double x[n];
	parse_real_data(n, x, arg); //Test the parsing itself.
	print_vector(n, x, "Parsed x");
	return assert_vector_equals(n, x, n_expc, x_expc);
}

/**
 * Automatically performs the parsing tests with the given real array "x", assuming there is no error.
 */
int test_parse_real_data_from_input_auto(int n, double* x) {
	char arg[20*(n+2)];
	int i, l = 0;
	for (i = 0; i < n; i++) {
		l += sprintf(arg+l, "%.16g ", x[i]);
	}
	return test_parse_real_data_from_input(arg, n, x);
}

/**
 * Make all tests on the parse_real_data() function.
 */
int test_parse_real_data() {
	print_test("parse_real_data()");
	double xa1[] = {0, 1, 2.1548, -3.5700, 4.12e12, -0.00001159, 1.5279845e5, -12.524e-15};
	double xa2[] = {124825487259, -129547852, -INFINITY, 1.e100, INFINITY};
	char arg1[] = "  1.1    2.2  \t 3.3 ";
	double x1_expc[] = {1.1, 2.2, 3.3};
	char arg2[] = "  0    inf  1.e-19  12458542478.";
	double x2_expc[] = {0, INFINITY, 1.e-19, 12458542478.};
	
	int nf = 0; //Number of failed tests.
	nf += test_parse_real_data_from_input_auto(sizeof(xa1)/sizeof(xa1[0]), xa1);
	nf += test_parse_real_data_from_input_auto(sizeof(xa2)/sizeof(xa2[0]), xa2);
	nf += test_parse_real_data_from_input(arg1, sizeof(x1_expc)/sizeof(x1_expc[0]), x1_expc);
	nf += test_parse_real_data_from_input(arg2, sizeof(x2_expc)/sizeof(x2_expc[0]), x2_expc);
	if (nf == 0) {
		printf(STR_PASS "parse_real_data(): All tests passed.\n");
	}
	return nf;
}

/**************************************
 * TESTING COMPLEX VECTORS UTILITIES:
 *************************************/
/**
 * Test the cnorm() function against some expected results.
 */
int test_cnorm_from_input(int n, dcomplex* z, double normz_expc) {
	print_cvector(n, z, "Input z");
	double normz = cnorm(n, z);
	int failed = (normz != normz_expc);
	if (failed) {
		printf(STR_FAIL "Invalid norm=%g, expected=%g \n", normz, normz_expc);
	}
	else {
		printf(STR_PASS "Norm(z) = %g \n", normz);
	}
	return failed;
}
/**
 * Test the cnorm() function.
 */
int test_cnorm() {
	print_test("cnorm()");
	dcomplex z1[] = {1, 2, 3, 4};
	double normz1_expc = sqrt(30.);
	dcomplex z2[] = {1+I, 2+2*I, 3+3*I, 4+4*I};
	double normz2_expc = sqrt(60.);
	dcomplex z3[] = {0, 0, 0, 0, 0};
	double normz3_expc = 0;
	
	int nf = 0; //Number of failed tests.
	nf += test_cnorm_from_input(sizeof(z1)/sizeof(z1[0]), z1, normz1_expc);
	nf += test_cnorm_from_input(sizeof(z2)/sizeof(z2[0]), z2, normz2_expc);
	nf += test_cnorm_from_input(sizeof(z3)/sizeof(z3[0]), z3, normz3_expc);
	if (nf == 0) {
		printf(STR_PASS "cnorm(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the cnormalize() function.
 */
int test_cnormalize() {
	print_test("cnormalize()");
	dcomplex x[] = {1e-9-2e-9*I, 1-4*I, 3+I, -2-3*I, 5-I, -1.};
	int n = sizeof(x)/sizeof(x[0]);
	int i, failed = 0;
	double err;
	
	dcomplex y[n], fac[n], factor;
	copy_cvector(n, x, y);
	print_cvector(n, x, "Original x");
	printf("[INFO] Norm(x) = %g \n", cnorm(n, x));
	
	factor = cnormalize(n, y); //Normalize "y" but not "x".
	print_cvector(n, y, "Normalized x");
	double normy = cnorm(n, y);
	if (fabs(normy - 1.) > TOLER) {
		failed = 1;
		printf(STR_FAIL "Invalid norm: Norm(x) = %g \n", normy);
	}
	else {
		printf(STR_PASS "Correct norm: Norm(x) = %g \n", normy);
	}
	for (i = 0; i < n; i++) {//Check the scaling factor.
		fac[i] = x[i]/y[i];
		err = cabs(fac[i]/factor - 1.); //Relative error.
		if (err > TOLER) {
			failed = 1; //Declare the test as failed.
			printf(STR_FAIL "Factor fac[%d]=%g%+gi is too far from expected %g%+gi (err=%g > %g).\n",
				i, creal(fac[i]), cimag(fac[i]), creal(factor), cimag(factor), err, TOLER);
		}
	}
	print_cvector(n, fac, "Factors");
	printf("[INFO] Returned factor = %g%+gi \n", creal(factor), cimag(factor));
	if (failed == 0) {
		printf(STR_PASS "cnormalize(): All tests passed.\n");
	}
	return failed;
}

/**
 * Test the distance between two complex vectors.
 */
int test_cdistance() {
	print_test("cdistance()");
	int failed = 0;
	dcomplex z[] = {2+I, 0, -3-2*I, 6+3*I, -2-I, 4};
	int n = sizeof(z)/(2*sizeof(z[0])); //n=3.
	print_cvector(n, z,   "x");
	print_cvector(n, z+n, "y");
	double dist = cdistance(n, z, z+n);
	double dist_expc = sqrt(78.);
	printf("[INFO] Distance |x-y| = %g\n", dist);
	if (fabs(dist/dist_expc - 1.) > TOLER) {
		failed = 1;
		printf(STR_FAIL "Invalid distance=%g, expected=%g.\n", dist, dist_expc);
	}
	dist = cdistance(n, z, z);
	dist_expc = 0.;
	printf("[INFO] Distance |x-x| = %g\n", dist);
	if (dist > TOLER) {
		failed = 1;
		printf(STR_FAIL "Invalid distance=%g, expected=%g.\n", dist, dist_expc);
	}
	if (failed == 0) {
		printf(STR_PASS "cdistance(): All tests passed.\n");
	}
	return failed;
}

/**
 * Test the initialization of constant complex vector with unit norm.
 */
int test_constant_unit_cvector() {
	print_test("constant_unit_cvector()");
	int n = 5;
	dcomplex z[n];
	constant_unit_cvector(n, z);
	print_cvector(n, z, "Constant z");
	
	double normz = cnorm(n, z);
	int failed = (fabs(normz - 1.) > TOLER);
	if (failed) {
		printf(STR_FAIL "Invalid norm: Norm(z) = %g\n", normz);
	}
	else {
		printf(STR_PASS "Correct norm: Norm(z) = %g\n", normz);
	}
	return failed;
}

/**
 * Test the sorting of complex vectors.
 */
int test_sort_cvector_from_input(int n, dcomplex* z, double (*f)(dcomplex)) {
	print_cvector(n, z, "Input  z");
	sort_cvector(n, z, f);
	print_cvector(n, z, "Sorted z");
	int i, failed = 0;
	for (i = 1; i < n; i++) {
		if (f(z[i-1]) > f(z[i])) {//The array must be sorted in increasing order by convention.
			failed = 1;
			printf(STR_FAIL "Invalid order between f(z[%d])=%g and f(z[%d])=%g, expected f(z[i-1]) < f(z[i]).\n", i-1, f(z[i-1]), i, f(z[i]));
		}
	}
	return failed;
}

/**
 * Make all tests of the sort_cvector() function.
 */
int test_sort_cvector() {
	print_test("sort_cvector()");
	dcomplex z1[] = {6.35-1.54*I, -7.77-5.05*I, 5.79+9.54*I, -6.24+6.5*I, -5.17+8.51*I, -8.69+1.56*I,
		0.84-4.14*I, -5.38-5.84*I, -2.08+1.61*I, 4.01-7.42*I, -5.76-3.87*I, 4.97+4.24*I};
	double f1(dcomplex z) {//Defines the sort function.
		return creal(z);
	};
	dcomplex z2[] = {0, 0, 1-I, 2-I, 4, 3, 0, 1+I, 1};
	double f2(dcomplex z) {
		return cabs(z);
	};
	dcomplex z3[] = {4.44-4.12*I, -7.81-6.7*I, -0.59+2.03*I, 0.71+5.08*I, 1.66+5.42*I, 0.};
	double f3(dcomplex z) {
		return cimag(z);
	};
	int nf = 0;
	nf += test_sort_cvector_from_input(sizeof(z1)/sizeof(z1[0]), z1, f1);
	nf += test_sort_cvector_from_input(sizeof(z2)/sizeof(z2[0]), z2, f2);
	nf += test_sort_cvector_from_input(sizeof(z3)/sizeof(z3[0]), z3, f3);
	if (nf == 0) {
		printf(STR_PASS "sort_cvector(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the deletion of duplicate entries in a complex vector.
 */
int test_uniq_cvector_from_input(int n, dcomplex* z, double toler, int n_expc, dcomplex* z_expc) {
	//int n_bak = n; //Backup the original number of components.
	print_cvector(n, z, "Input z");
	uniq_cvector(&n, z, toler);
	print_cvector(n, z, "Uniqs z");
	//print_cvector(n_bak, z, "Whole z");
	return assert_cvector_equals(n, z, n_expc, z_expc);
}

/**
 * Make all the test on the uniq_cvector() function.
 */
int test_uniq_cvector() {
	print_test("uniq_cvector()");
	dcomplex z1[] = {-0.43-7.21*I, -0.43-7.21*I, -9.83-6.39*I, -3.06+0.57*I, -9.83-6.39*I, -3.06+0.57*I, -3.06+0.57*I, -0.43-7.21*I, -0.43-7.21*I, -9.83-6.39*I};
	dcomplex z1_expc[] = {-0.43-7.21*I, -9.83-6.39*I, -3.06+0.57*I};
	double toler1 = 1e-15; //Tolerate nearly zero relative difference.
	
	dcomplex z2[] = {-3.77+2.54*I, 1.32-7.41*I, -3.776+2.549*I, 1.314-7.401*I};
	dcomplex z2_expc[] = {-3.77+2.54*I, 1.32-7.41*I};
	double toler2 = 0.01; //Tolerate 1% relative difference.
	
	dcomplex z3[] = {0, 0, 3+I, 3+I, 0, 3+I, 2, 2, 2, 3.01+I, 3+I, 2, 0.05, -0.03, -0.02+0.01*I, 3+I, 2, 2.1, -0.01, 2, 2, 3+I};
	dcomplex z3_expc[] = {0, 3+I, 2};
	double toler3 = 0.1; //Tolerate 10% relative difference.
	
	int nf = 0;
	nf += test_uniq_cvector_from_input(sizeof(z1)/sizeof(z1[0]), z1, toler1, sizeof(z1_expc)/sizeof(z1_expc[0]), z1_expc);
	nf += test_uniq_cvector_from_input(sizeof(z2)/sizeof(z2[0]), z2, toler2, sizeof(z2_expc)/sizeof(z2_expc[0]), z2_expc);
	nf += test_uniq_cvector_from_input(sizeof(z3)/sizeof(z3[0]), z3, toler3, sizeof(z3_expc)/sizeof(z3_expc[0]), z3_expc);
	if (nf == 0) {
		printf(STR_PASS "uniq_cvector(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the filtering of the complex vector entries using the filter_cvector() function.
 */
int test_filter_cvector_from_input(int n, dcomplex* z, int* flag, int n_expc, dcomplex* z_expc) {
	int i;
	print_cvector(n, z, "Input z");
	printf("[INFO] Input flags = {");
	for (i = 0; i < n; i++) {//Prints the flag array.
		printf("%d", flag[i]);
		if (i != n-1)
			printf(", ");
	}
	printf("}.\n");
	//int n_bak = n; //Backup the original number of components.
	filter_cvector(&n, z, flag);
	print_cvector(n, z, "Result z");
	//print_cvector(n_bak, z, "Whole  z");
	return assert_cvector_equals(n, z, n_expc, z_expc);
}

/**
 * Make all the tests for the filter_cvector() function.
 */
int test_filter_cvector() {
	print_test("filter_cvector()");
	int flag1[] = {0, 1, 1, -1, 1, 0, 1, 2, 0, 1};
	dcomplex z1[] = {6.6+8.6*I, 4+9*I, 2.4+2.9*I, 3-3.8*I, 8.3+2.2*I, 4.4-8.6*I, 7.5+9.5*I, 9.9+0.9*I, -8.3+2.8*I, -0.8+1.7*I};
	dcomplex z1_expc[] = {4+9*I, 2.4+2.9*I, 8.3+2.2*I, 7.5+9.5*I, -0.8+1.7*I};
	
	int flag2[] = {0, 1, 0, 0, 1, 0};
	dcomplex z2[] = {1, 2, 3, 4, 5, 6};
	dcomplex z2_expc[] = {2, 5};
	
	int flag3[] = {-1, -1, -1};
	dcomplex z3[] = {4+9*I, 4.4-8.6*I, -8.3+2.8*I};
	dcomplex z3_expc[] = {};
	
	int nf = 0; //Number of failed tests.
	nf += test_filter_cvector_from_input(sizeof(z1)/sizeof(z1[0]), z1, flag1, sizeof(z1_expc)/sizeof(z1_expc[0]), z1_expc);
	nf += test_filter_cvector_from_input(sizeof(z2)/sizeof(z2[0]), z2, flag2, sizeof(z2_expc)/sizeof(z2_expc[0]), z2_expc);
	nf += test_filter_cvector_from_input(sizeof(z3)/sizeof(z3[0]), z3, flag3, sizeof(z3_expc)/sizeof(z3_expc[0]), z3_expc);
	if (nf == 0) {
		printf(STR_PASS "filter_cvector(): All tests passed.\n");
	}
	return nf;
}

/**
 * Test the complex vector parsing function.
 */
int test_parse_complex_data_from_input(char* arg, int n_expc, dcomplex* z_expc) {
	printf("[INFO] Input string = '%s'\n", arg);
	int n = count_complex_data(arg);  //Test the estimate of arguments.
	dcomplex z[n];
	parse_complex_data(n, z, arg); //Test the parsing itself.
	print_cvector(n, z, "Parsed z");
	return assert_cvector_equals(n, z, n_expc, z_expc);
}

/**
 * Automatically performs the parsing tests with the given complex array "z", assuming there is no error.
 */
int test_parse_complex_data_from_input_auto(int n, dcomplex* z) {
	char arg[40*(n+2)];
	int i, l = 0;
	for (i = 0; i < n; i++) {
		l += sprintf(arg+l, "%.16g%+.16gi ", creal(z[i]), cimag(z[i]));
	}
	return test_parse_complex_data_from_input(arg, n, z);
}

/**
 * Make all tests on the parse_complex_data() function.
 */
int test_parse_complex_data() {
	print_test("parse_complex_data()");
	dcomplex za1[] = {1.3+2.5*I, 3+I, -1.8-4.7*I, 0, 3.34e-7-1.15e-9*I, 10.3+12.5*I};
	dcomplex za2[] = {1, 2, 1215548759+8524579624*I, 11*I, INFINITY, 10.3+12.5*I};
	char arg1[] = "1.3+2.5i  \t  3.34e-7+1.15e-9i  0+11i -3.9-9.8i\t lol";
	dcomplex z1_expc[] = {1.3+2.5*I, 3.34e-7+1.15e-9*I, 11*I, -3.9-9.8*I};
	char arg2[] = "\t1+0i  0+0i   2.6-5.4i     3+1i";
	dcomplex z2_expc[] = {1, 0, 2.6-5.4*I, 3+I};
	
	int nf = 0; //Number of failed tests.
	nf += test_parse_complex_data_from_input_auto(sizeof(za1)/sizeof(za1[0]), za1); //Automatic tests.
	nf += test_parse_complex_data_from_input_auto(sizeof(za2)/sizeof(za2[0]), za2);
	nf += test_parse_complex_data_from_input(arg1, sizeof(z1_expc)/sizeof(z1_expc[0]), z1_expc);
	nf += test_parse_complex_data_from_input(arg2, sizeof(z2_expc)/sizeof(z2_expc[0]), z2_expc);
	if (nf == 0) {
		printf(STR_PASS "parse_complex_data(): All tests passed.\n");
	}
	return nf;
}

/**
 * Main function.
 */
int main(int argc, char** argv) {
	
	int nf = 0; //Current number of failed tests.
	
	/*** Real vector tests ***/
	nf += test_norm();
	nf += test_distance();
	nf += test_sort_vector();
	nf += test_uniq_vector();
	nf += test_quantile();
	nf += test_strtod();
	nf += test_parse_real_data();
	
	/*** Complex vector tests ***/
	nf += test_cnorm();
	nf += test_cnormalize();
	nf += test_cdistance();
	nf += test_constant_unit_cvector();
	nf += test_sort_cvector();
	nf += test_uniq_cvector();
	nf += test_filter_cvector();
	nf += test_parse_complex_data();
	
	printf("====== TEST RESULTS ======\n");
	if (nf == 0) {
		printf(STR_PASS "All tests are successful.\n");
	}
	else {
		printf(STR_FAIL "There are %d failed tests.\n", nf);
	}
	return 0;
}
