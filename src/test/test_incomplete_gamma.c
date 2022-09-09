/****
 * @date Created on 2021-08-24 at 13:34:40 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to systematically test the computations of the incomplete gamma functions.
 ***/
#include <stdio.h>    /* Standard Library of Input and Output */
#include <math.h>     /* Standard Library for Mathematical Functions */
#include "incomplete_gamma.h" /* Import the Incomplete Gamma and Related Functions */
#include "tag_strings.h"
#include "assertion.h"

#define TOLER  1.e-13  //General tolerance on results.


/**
 * Test the confluent hypergeometric function 1F1(a,b,z) computed from its series representation for some values.
 */
int test_lower_gamma() {
	print_test("gamma_lower()");
	double lowgamma_data[][4] = {//Use Mathematica or any other reliable software to generate this data set.
		//{a, z, f},
		{1.5, 0.1, 0.019860967741930694778},
		{1.0,   2., 0.86466471676338730811},
		{1.5,   3., 0.78731493881798064404},
		{2.0,   5., 0.95957231800548719742},
		{2.5,  10., 1.3276790708673575604},
		{1.5,  20., 0.88622691600993006630},
		{1.5,  50., 0.88622692545275801365},
		{1.5, 150., 0.88622692545275801365},
		{5.5,   2., 1.5746265342113649474},
		{5.5,  10., 49.969520938049710198},
		{5.5,  20., 52.340905215160154393},
		{5.5,  50., 52.342777784553510834},
		{5.5, 150., 52.342777784553520181},
		{10.,  20., 361067.26478156134277},
		{10., 150., 362880.00000000000000}
	};
	int ntest = sizeof(lowgamma_data)/sizeof(lowgamma_data[0]);
	printf("[INFO] Number of tests: %d\n", ntest);
	double a, z, f, f_expc, err;
	char* tag;
	int i, nfail = 0;  //Number of failed tests.
	for (i = 0; i < ntest; i++) {//Loop on the tests.
		a = lowgamma_data[i][0];
		z = lowgamma_data[i][1];
		f = lower_gamma(a, z);
		f_expc = lowgamma_data[i][2];
		err = fabs(f/f_expc - 1.);
		if (err > TOLER) {
			tag = STR_FAIL;
			nfail++;
		}
		else {
			tag = STR_PASS;
		}
		printf("%s#%2d | a=%6g, z=%6g, g=%12g (expected g=%12g, err=%g) \n", tag, i, a, z, f, f_expc, err);
	}
	return nfail;
}

/**
 * Test the confluent hypergeometric function 1F1(a,b,z) computed from its series representation for some values.
 */
int test_hypergeom_1f1_series() {
	print_test("hypergeom_1f1_series()");
	printf(STR_WARN "This is an internal function which may not work everywhere...\n");
	double hyp1f1_data[][4] = {//Use Mathematica or any reliable software to generate this data set.
		//{a, b, z, f},
		{1., 1.0,   1., 2.71828182845904524},
		{1., 1.5,   3., 10.1300112010004876},
		{1., 2.0,   5., 29.4826318205153207},
		{1., 2.5,  -1., 0.692880739630847371},
		{1., 3.0,  15., 29057.7899775298723},
		{2., 3.5,   2., 3.47799357209040213},
		{3., 5.0,  15., 132504.162297536218},
		{1., 4.0, 400., 4.8951278341538849537e166},
	};
	int ntest = sizeof(hyp1f1_data)/sizeof(hyp1f1_data[0]);
	printf("[INFO] Number of tests: %d\n", ntest);
	double a, b, z, f, f_expc, err;
	char* tag;
	int i, nfail = 0;  //Number of failed tests.
	for (i = 0; i < ntest; i++) {//Loop on the tests.
		a = hyp1f1_data[i][0];
		b = hyp1f1_data[i][1];
		z = hyp1f1_data[i][2];
		f = hypergeom_1f1_series(a, b, z);
		f_expc = hyp1f1_data[i][3];
		err = fabs(f/f_expc - 1.);
		tag = (err > TOLER)? STR_WARN : STR_PASS;
		printf("%s#%2d | a=%6g, b=%6g, z=%6g, 1f1=%12g (expected 1f1=%12g, err=%g) \n", tag, i, a, b, z, f, f_expc, err);
	}
	return nfail;
}

/**
 * Main function
 */
int main(int argc, char** argv) {
	
	int nf = 0; //Current number of failed tests.
	
	nf += test_lower_gamma();
	nf += test_hypergeom_1f1_series();
	
	printf("====== TEST RESULTS ======\n");
	if (nf == 0) {
		printf(STR_PASS "All tests are successful.\n");
	}
	else {
		printf(STR_FAIL "There are %d failed tests.\n", nf);
	}
	return 0;
}
