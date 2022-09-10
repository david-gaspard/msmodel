/****
 * @date Created on 2021-04-29 at 17:02:16 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to systematically test the computations of Bessel functions and Green functions.
 ***/
#include <stdio.h>    /* Standard Library of Input and Output */
#include <math.h>     /* Standard Library for Mathematical Functions */
#include "tag_strings.h"
#include "assertion.h"
#include "green_bessel.h"
#include "test_green_bessel_data.h" /* Import the data to the the Bessel and Green functions */

#define TOLER  1.e-12  //Global tolerance threshold on the resulting values used for test.

/**
 * Displays the volume and surface area of the unit d-ball for the first dimensions.
 */
void show_uball_volsurf(int dmax) {
	int d;
	printf("[INFO] Volume  = {");
	for (d = 0; d <= dmax; d++) {
		printf("%d: %g", d, uball_volume(d));
		if (d != dmax) {
			printf(",\t");
		}
	}
	printf("},\n[INFO] Surface = {");
	for (d = 0; d <= dmax; d++) {
		printf("%d: %g", d, uball_surface(d));
		if (d != dmax) {
			printf(",\t");
		}
	}
	printf("}.\n");
}

/**
 * Returns the relative error on the computed K_nu(z) from arguments on the format [nu, Re(z), Im(z), Re(K_nu(z)), Im(K_nu(z))].
 */
double bessel_k_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "nu".
	if (fracp != 0 && fracp != 0.5 && fracp != -0.5) {
		printf("[WARN] nu=%g, only integer and half-integer orders are supported for K_nu(z).\n", arg[0]);
	}
	int twonu = (int)floor(2*arg[0]); //Computes "twonu" from "nu".
	dcomplex z = arg[1] + I*arg[2];
	dcomplex f = bessel_k(twonu, z);
	dcomplex expc_f = arg[3] + I*arg[4];
	double err = cabs(f/expc_f - 1.);
	if (err > TOLER) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] nu=%g, z=%5g%+5gi, K_nu(z)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			twonu/2., creal(z), cimag(z), creal(f), cimag(f), creal(expc_f), cimag(expc_f), err);
	}
	return err;
}

/**
 * Returns the relative error on the computed G(d,k,r) from arguments on the format [d, Re(k), Im(k), Re(G(d,k,r)), Im(G(d,k,r))].
 */
double free_green_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "d".
	if (fracp != 0) {
		printf("[WARN] d=%g, only integer dimensions are supported for G(d,k,r).\n", arg[0]);
	}
	int d = (int)floor(arg[0]);
	dcomplex k = arg[1] + I*arg[2];
	double r = 1.;  //All the test assume r=1.
	dcomplex g = free_green(d, k, r);
	dcomplex g_expc = arg[3] + I*arg[4];
	double err = cabs(g/g_expc - 1.);
	if (err > TOLER) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] d=%d, k=%5g%+5gi, G(d,k,r)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			d, creal(k), cimag(k), creal(g), cimag(g), creal(g_expc), cimag(g_expc), err);
	}
	return err;
}

/**
 * Returns the relative error on the computed I(d,k,r) from arguments on the format [d, Re(k), Im(k), Re(I(d,k,r)), Im(I(d,k,r))].
 */
double free_green_imag_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "d".
	if (fracp != 0) {
		printf("[WARN] d=%g, only integer dimensions are supported for I(d,k,r).\n", arg[0]);
	}
	int d = (int)floor(arg[0]);
	dcomplex k = arg[1] + I*arg[2];
	double r = 1.;  //All the test assume r=1.
	dcomplex i = free_green_imag(d, k, r);
	dcomplex i_expc = arg[3] + I*arg[4];
	double err = cabs(i/i_expc - 1.);
	if (err > TOLER) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] d=%d, k=%5g%+5gi, I(d,k,r)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			d, creal(k), cimag(k), creal(i), cimag(i), creal(i_expc), cimag(i_expc), err);
	}
	return err;
}

/**
 * Returns the relative error on the computed dG/dk from arguments on the format [d, Re(k), Im(k), Re(dG/dk), Im(dG/dk)].
 */
double dk_free_green_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "d".
	if (fracp != 0) {
		printf("[WARN] d=%g, only integer dimensions are supported for dG/dk(d,k,r).\n", arg[0]);
	}
	int d = (int)floor(arg[0]);
	dcomplex k = arg[1] + I*arg[2];
	double r = 1.;  //All the test assume r=1.
	dcomplex dg = dk_free_green(d, k, r);
	dcomplex dg_expc = arg[3] + I*arg[4];
	double err = cabs(dg/dg_expc - 1.);
	if (err > TOLER) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] d=%d, k=%5g%+5gi, dG/dk(d,k,r)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			d, creal(k), cimag(k), creal(dg), cimag(dg), creal(dg_expc), cimag(dg_expc), err);
	}
	return err;
}

/**
 * Returns the relative error on the ratio of Bessel functions J_(nu+1)(z)/J_nu(z) from arguments on the format [nu, Re(z), Im(z), Re(J_nu+1/J_nu), Im(J_nu+1/J_nu)].
 */
double bessel_j_ratio_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "nu".
	if (fracp != 0 && fracp != 0.5 && fracp != -0.5) {
		printf("[WARN] nu=%g, only integer and half-integer orders are supported for J_nu(z).\n", arg[0]);
	}
	int twonu = (int)floor(2*arg[0]); //Computes "twonu" from "nu".
	dcomplex z = arg[1] + I*arg[2];
	dcomplex f = bessel_j_ratio(twonu, z);
	dcomplex expc_f = arg[3] + I*arg[4];
	double err = cabs(f/expc_f - 1.);
	if (err > TOLER) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] nu=%g, z=%5g%+5gi, J_(nu+1)(z)/J_nu(z)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			twonu/2., creal(z), cimag(z), creal(f), cimag(f), creal(expc_f), cimag(expc_f), err);
	}
	return err;
}

/**
 * Returns the relative error on the ratio of Bessel functions K_(nu+1)(z)/K_nu(z) from arguments on the format [nu, Re(z), Im(z), Re(K_nu+1/K_nu), Im(K_nu+1/K_nu)].
 */
double bessel_k_ratio_error(double* arg) {
	double fracp = fmod(arg[0], 1.);  //Fractional part of "nu".
	if (fracp != 0 && fracp != 0.5 && fracp != -0.5) {
		printf("[WARN] nu=%g, only integer and half-integer orders are supported for K_nu(z).\n", arg[0]);
	}
	int twonu = (int)floor(2*arg[0]); //Computes "twonu" from "nu".
	dcomplex z = arg[1] + I*arg[2];
	dcomplex f = bessel_k_ratio(twonu, z);
	dcomplex expc_f = arg[3] + I*arg[4];
	double err = cabs(f/expc_f - 1.);
	if (err > TOLER || err != err) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] nu=%g, z=%5g%+5gi, K_(nu+1)(z)/K_nu(z)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			twonu/2., creal(z), cimag(z), creal(f), cimag(f), creal(expc_f), cimag(expc_f), err);
	}
	return err;
}

/**
 * Returns the relative error on the regularized hypergoemetric function 0F1(a,z) from arguments on the format [a, Re(z), Im(z), Re(0F1(a,z)), Im(0F1(a,z))].
 */
double hypergeom_0f1_reg_series_error(double* arg) {
	double a = arg[0];
	dcomplex z = arg[1] + I*arg[2];
	dcomplex f = hypergeom_0f1_reg_series(a, z);
	dcomplex expc_f = arg[3] + I*arg[4];
	double err = cabs(f/expc_f - 1.);
	if (err > TOLER || err != err) {//Only shows the errors above the tolerance threshold.
		printf("[WARN] a=%g, z=%5g%+5gi, 0F1(a,z)=%10g%+10gi (expected=%10g%+10gi), err=%g \n",
			a, creal(z), cimag(z), creal(f), cimag(f), creal(expc_f), cimag(expc_f), err);
	}
	return err;
}

/**
 * Test any function "f(z)" against the expected results from the data array "test".
 * The relative errors on "f(z)" are given by the function "funcerr()" which takes an array of arguments.
 */
void test_function(const char* funcname, double (*funcerr)(double*), int ntest, double (*test)[5]) {
	print_test(funcname);
	printf("[INFO] Number of tests: %d\n", ntest);
	double err, sumerr = 0., logerr = 0.;
	int i;
	for (i = 0; i < ntest; i++) {//Loop on the tests.
		err = funcerr(test[i]);
		sumerr += err;
		logerr += (err != 0.) ? log(err) : 0.;
	}
	printf("[INFO] Relative error (AM) = %g.\n", sumerr/ntest);
	printf("[INFO] Relative error (GM) = %g.\n", exp(logerr/ntest));
}

/**
 * Main function
 */
int main(int argc, char** argv) {
	
	show_uball_volsurf(7);
	
	int ntest;
	
	ntest = sizeof(BESSEL_K_DATA)/sizeof(BESSEL_K_DATA[0]);
	test_function("bessel_k()", bessel_k_error, ntest, BESSEL_K_DATA);
	
	ntest = sizeof(GREEN_DATA)/sizeof(GREEN_DATA[0]);
	test_function("free_green()", free_green_error, ntest, GREEN_DATA);
	
	ntest = sizeof(HYP0F1_SERIES_DATA)/sizeof(HYP0F1_SERIES_DATA[0]);
	test_function("hypergeom_0f1_reg_series()", hypergeom_0f1_reg_series_error, ntest, HYP0F1_SERIES_DATA);
	
	ntest = sizeof(GREEN_IMAG_DATA)/sizeof(GREEN_IMAG_DATA[0]);
	test_function("free_green_imag()", free_green_imag_error, ntest, GREEN_IMAG_DATA);
	
	ntest = sizeof(DK_GREEN_DATA)/sizeof(DK_GREEN_DATA[0]);
	test_function("dk_free_green()", dk_free_green_error, ntest, DK_GREEN_DATA);
	
	ntest = sizeof(BESSEL_J_RATIO_DATA)/sizeof(BESSEL_J_RATIO_DATA[0]);
	test_function("bessel_j_ratio()", bessel_j_ratio_error, ntest, BESSEL_J_RATIO_DATA);
	
	ntest = sizeof(BESSEL_K_RATIO_DATA)/sizeof(BESSEL_K_RATIO_DATA[0]);
	test_function("bessel_k_ratio()", bessel_k_ratio_error, ntest, BESSEL_K_RATIO_DATA);
	
	return 0;
}
