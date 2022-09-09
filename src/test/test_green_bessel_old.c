/**************************************************************************************************
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created at 2020-07-25 at 18:13 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Test the output of the Bessel and Green functions
 **************************************************************************************************/
#include <stdio.h>    /* Standard Library of Input and Output */
#include <math.h>     /* Standard Library for Mathematical Functions */
#include <time.h>     /* Standard Library of CPU Time Measurement */
#include "green_bessel.h"
/**
 * Define some preprocessor constants:
 */
#define NTEST   1000  //Number of test while computing the taken time. At least 10^5 to get enough statistics.
/******************************************
 * TESTING FUNCTIONS
 ******************************************/
/**
 * Speed test of the Bessel K_nu(z) function.
 */
void bessel_k_int_speed_test(int nu, dcomplex z, FILE* fp) {
	dcomplex kres;
	int i;
	clock_t t = clock();
	for (i = 0; i < NTEST; i++) {
		kres = bessel_k_int(nu, z);
	}
	double dt = ((double)(clock() - t)/CLOCKS_PER_SEC)/NTEST; //CPU time per function in seconds.
	printf("[INFO] K_%d(%f %+fi) =\t %.16e %+.16ei \t[avgdt=%es, ntest=%d]\n", nu, creal(z), cimag(z), creal(kres), cimag(kres), dt, NTEST);
	fprintf(fp, "%e\t%e\t%.16e\t%.16e\t%e\n", creal(z), cimag(z), creal(kres), cimag(kres), dt);
}

/**
 * Checks the implementation of the Bessel K_nu(z) function for integer nu, and test the performances.
 * The function scans the results of K_nu(z) at different key points in the z-plane, and outputs the results to a file.
 */
void bessel_k_int_full_test(int nu) {
	dcomplex z;
	char filename[50];
	sprintf(filename, "test_besselk_int_nu%d.tab", nu);
	FILE* fp = fopen(filename, "w");
	int i;
	for (i = 0; i < 10; i++) {
		z = 2.0*i+1.0;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 1; i < 10; i++) {
		z = 50.0*i;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 0; i < 10; i++) {
		z = (2.0*i+1.0)*I;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 1; i < 10; i++) {
		z = (50.0*i)*I;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 0; i < 10; i++) {
		z = -(2.0*i+1.0) - 0.1*I;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 1; i < 10; i++) {
		z = -50.0*i - 1.0*I;
		bessel_k_int_speed_test(nu, z, fp);
	}
	for (i = 0; i < 10; i++) {
		z = (2.0*i+1.0)*(-1.0-1.0*I);
		bessel_k_int_speed_test(nu, z, fp);
	}
	fclose(fp);
}

/**
 * Speed test of the free Green function G(d,k,r).
 */
void free_green_speed_test(int d, dcomplex k, double r, FILE* fp) {
	dcomplex g;
	int i;
	clock_t t = clock();
	for (i = 0; i < NTEST; i++) {
		g = free_green(d, k, r);
	}
	double dt = ((double)(clock() - t)/CLOCKS_PER_SEC)/NTEST; //CPU time per function in seconds.
	printf("[INFO] G%d(k=%f %+fi, r=%f) =\t %.16e %+.16ei \t[avgdt=%es, ntest=%d]\n", d, creal(k), cimag(k), r, creal(g), cimag(g), dt, NTEST);
	fprintf(fp, "%e\t%e\t%.16e\t%.16e\t%e\n", creal(k), cimag(k), creal(g), cimag(g), dt);
}

/**
 * Checks the implementation of the Green function G(d,k,r) function for integer d, and test the performances.
 * The function scans the results of G(d,k,r) at different key points in the k-plane, assuming r=1, and outputs the results to a file.
 */
void free_green_full_test(int d) {
	dcomplex k;
	double r = 1.0;
	char filename[50];
	sprintf(filename, "test_green_d%d.tab", d);
	FILE* fp = fopen(filename, "w");
	int i;
	for (i = 0; i < 10; i++) {
		k = 2.0*i+1.0;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = 50.0*i;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = (2.0*i+1.0)*I;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = (50.0*i)*I;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = -(2.0*i+1.0) - 0.1*I;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = -50.0*i - 1.0*I;
		free_green_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = (2.0*i+1.0)*(-1.0-1.0*I);
		free_green_speed_test(d, k, r, fp);
	}
	fclose(fp);
}


/**
 * Speed test of the imaginary part of the free Green function, -Im[G_k(r)].
 */
void free_green_imag_speed_test(int d, dcomplex k, double r, FILE* fp) {
	dcomplex q;
	int i;
	clock_t t = clock();
	for (i = 0; i < NTEST; i++) {
		q = free_green_imag(d, k, r);
	}
	double dt = ((double)(clock() - t)/CLOCKS_PER_SEC)/NTEST; //CPU time per function in seconds.
	printf("[INFO] Q%d(k=%f %+fi, r=%f) =\t %.16e %+.16ei \t[avgdt=%es, ntest=%d]\n", d, creal(k), cimag(k), r, creal(q), cimag(q), dt, NTEST);
	fprintf(fp, "%e\t%e\t%.16e\t%.16e\t%e\n", creal(k), cimag(k), creal(q), cimag(q), dt);
}

/**
 * Checks the implementation of the imaginary part of the Green function, -Im[G_k(r)], for integer d, and test the performances.
 * The function scans the results of free_green_imag() at different key points in the k-plane, assuming r=1, and outputs the results to a file.
 */
void free_green_imag_full_test(int d) {
	dcomplex k;
	double r = 1.0;
	char filename[50];
	sprintf(filename, "test_green_imag_d%d.tab", d);
	FILE* fp = fopen(filename, "w");
	int i;
	for (i = 0; i < 10; i++) {
		k = 2.0*i+1.0;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = 50.0*i;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = (2.0*i+1.0)*I;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = (50.0*i)*I;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = -(2.0*i+1.0) - 0.1*I;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 1; i < 10; i++) {
		k = -50.0*i - 1.0*I;
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = (2.0*i+1.0)*(1.0-1.0*I);
		free_green_imag_speed_test(d, k, r, fp);
	}
	for (i = 0; i < 10; i++) {
		k = (2.0*i+1.0)*(-1.0-1.0*I);
		free_green_imag_speed_test(d, k, r, fp);
	}
	fclose(fp);
}

/**
 * Main function used for testing the implementation and the performance of the Bessel function computation.
 */
int default_main(int argc, char** argv) {
	
	int d;
	
	if (argc > 1) {
		int value;
		if (sscanf(argv[1], "%d", &value) == 1) {
			d = value;
		}
		else {
			printf("[ERROR] Incorrect dimension '%s', expecting a relatively small integer such as 1, 2, 3.\n", argv[1]);
			return 1;
		}
	}
	else {
		printf("[ERROR] No dimension found, aborting...\n[USAGE] %s d, where d is the integer dimension.\n", argv[0]);
		return 1;
	}
	if (d < 1) {
		printf("[ERROR] The dimension '%d' should be a strictly positive integer, aborting...\n", d);
		return 1;
	}
	
	//printf("[INFO] Testing the Green function G%d(k, r=1)...\n", d);
	//free_green_full_test(d);
	
	printf("[INFO] Testing the imaginary part of the Green function Q%d(k, r=1)...\n", d);
	free_green_imag_full_test(d);
	
	if (d%2 == 0) {
		int nu = (d - 2)/2; //Convert to the Bessel order.
		printf("[INFO] Even dimension detected, testing the Bessel K_%d(z)...\n", nu);
		bessel_k_int_full_test(nu);
	}
	return 0;
}

/**
 * Checks the implementation of the Bessel K_nu(z) function for integer or half-integer "nu", without testing the performances.
 */
void test_bessel_k() {
	char filename[] = "test_besselk.tab";
	FILE* fp = fopen(filename, "w");
	dcomplex z, k;
	double fracp; //Fractional part.
	int twonu;  //Stores the value of 2*nu.
	//Set of test points, all the combinations must be checked:
	double nu[] = {0, 0.5, 1, 1.5, 2, 3.5, 5., 10., 10.5, 20., 20.5, 50., 50.5};
	dcomplex zre[] = {-150., -10., -1., 0., 0.01, 1., 10., 150.};
	dcomplex zim[] = {-100., -10., -0.1, 0., 0.01, 1., 10., 100.};
	int inu, nulen = sizeof(nu)/sizeof(double);
	int izr, zrelen = sizeof(zre)/sizeof(dcomplex);
	int izi, zimlen = sizeof(zim)/sizeof(dcomplex);
	for (inu = 0; inu < nulen; inu++) {
		fracp = fmod(nu[inu], 1.);  //Fractional part of "nu".
		if (fracp != 0 && fracp != 0.5) {
			printf("[ERROR] Only integer and half-integer orders are supported, aborting...\n");
			return;
		}
		twonu = (int)(2*nu[inu]);  //Safe cast to integer 2*nu.
		for (izr = 0; izr < zrelen; izr++) {
			for (izi = 0; izi < zimlen; izi++) {
				z = zre[izr] + I*zim[izi];
				k = bessel_k(twonu, z);
				fprintf(fp, "%g\t%.16g\t%.16g\t%.17g\t%.17g\n", nu[inu], creal(z), cimag(z), creal(k), cimag(k));
			}
		}
	}
	fclose(fp);
}

/**
 * Test the computation of the ratio of Bessel functions, J_(nu+1)(z)/J_nu(z), and K_(nu+1)(z)/K_nu(z).
 */
void test_bessel_ratio() {
	FILE* fj = fopen("test_besselj_ratio.tab", "w");
	FILE* fk = fopen("test_besselk_ratio.tab", "w");
	dcomplex z, bjr, bkr;
	double fracp; //Fractional part.
	int twonu;  //Stores the value of 2*nu.
	//Set of test points, all the combinations must be checked:
	double nu[] = {0, 0.5, 1, 1.5, 2, 3.5, 5., 10., 10.5, 20., 20.5, 50., 50.5};
	dcomplex zre[] = {-150., -20., -10., -1., 0., 0.01, 1., 10., 20., 150.};
	dcomplex zim[] = {-100., -20., -10., -0.1, 0., 0.01, 1., 10., 20., 100.};
	int inu, nulen = sizeof(nu)/sizeof(double);
	int izr, zrelen = sizeof(zre)/sizeof(dcomplex);
	int izi, zimlen = sizeof(zim)/sizeof(dcomplex);
	for (inu = 0; inu < nulen; inu++) {
		fracp = fmod(nu[inu], 1.);  //Fractional part of "nu".
		if (fracp != 0 && fracp != 0.5) {
			printf("[ERROR] Only integer and half-integer orders are supported, aborting...\n");
			return;
		}
		twonu = (int)(2*nu[inu]);  //Safe cast to integer 2*nu.
		for (izr = 0; izr < zrelen; izr++) {
			for (izi = 0; izi < zimlen; izi++) {
				z = zre[izr] + I*zim[izi];
				bjr = bessel_j_ratio(twonu, z);
				bkr = bessel_k_ratio(twonu, z);
				fprintf(fj, "%g\t%.16g\t%.16g\t%.17g\t%.17g\n", nu[inu], creal(z), cimag(z), creal(bjr), cimag(bjr));
				fprintf(fk, "%g\t%.16g\t%.16g\t%.17g\t%.17g\n", nu[inu], creal(z), cimag(z), creal(bkr), cimag(bkr));
			}
		}
	}
	fclose(fj);
	fclose(fk);
}


/**
 * Test the computation of the imaginary part of the Green function, I(d,k,r) = -Im[G(d,k,r)].
 */
void test_green_imag() {
	FILE* fp = fopen("test_green_imag.tab", "w");
	dcomplex k, ig;
	//Set of test points, all the combinations must be checked:
	double r = 1.;
	int d[] = {1, 2, 3, 4, 5, 6, 7};
	dcomplex kre[] = {-100., -20., -10., -1., -0.1, 0., 0.01, 1., 10., 20., 100.};
	dcomplex kim[] = {-100., -20., -10., -1., -0.1, 0., 0.01, 1., 10., 20., 100.};
	int id, dlen = sizeof(d)/sizeof(int);
	int ikr, krelen = sizeof(kre)/sizeof(dcomplex);
	int iki, kimlen = sizeof(kim)/sizeof(dcomplex);
	for (id = 0; id < dlen; id++) {
		for (ikr = 0; ikr < krelen; ikr++) {
			for (iki = 0; iki < kimlen; iki++) {
				k = kre[ikr] + I*kim[iki];
				ig = free_green_imag(d[id], k, r);
				fprintf(fp, "%d\t%.16g\t%.16g\t%.17g\t%.17g\n", d[id], creal(k), cimag(k), creal(ig), cimag(ig));
			}
		}
	}
	fclose(fp);
}

/**
 * Main function.
 */
int main(int argc, char** argv) {
	
	test_bessel_k();
	test_bessel_ratio();
	test_green_imag();
	
	//return default_main(argc, argv);
	return 0;
}
