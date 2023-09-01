/****
 * @date Created on 2021-02-09 at 16:17:15 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to test the building and computation of the multi-scattering matrix.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tag_strings.h"
#include "assertion.h"
#include "complex_vector_util.h"
#include "medium_util.h"
#include "ms_matrix_util.h"
/**
 * Defines some macros:
 */
#define PI      3.1415926535897932385
#define TWOPI   6.2831853071795864769
/**
 * LAPACK's routine to solve a complex linear system using the LU decomposition method:
 */
void zgesv_(int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv, dcomplex* b, int* ldb, int* info);

/**
 * Test both the tracelog() function computing the trace of the logarithm of the matrix passed in argument,
 * and the build_matrix() function computing the multi-scattering matrix for some prescribed medium.
 */
void test_tracelog() {
	print_test("tracelog_lu()");
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {
		"dimension=2",
		"natom=50",
		"model=hardsphere 0.4",
		"shape=ball",
		"ratio=1."
	};
	uint64_t seed = 1; //Medium seed.
	
	parse_medium(med, sizeof(args)/sizeof(args[0]), args);
	print_param_medium(med);
	fill_medium(med, seed);
	//save_medium_points(med); exit(EXIT_SUCCESS); //Optionally save the atomic positions to test with an external program.
	
	//Test several values of "k" against expected values of log(det(M))/N computed in Wolfram Mathematica with the above settings for seed=1.
	dcomplex ks[] = {1. - 1.*I, 3. - 2.*I, 6. - 6.*I, 10. - 0.5*I, 10. - 10.*I};
	dcomplex logdet_expc[] = {-0.7384505297011151 - 0.01552793781658203*I,
		3.198589101932406 - 0.04336958744130083*I,
		22.04929061014481851 - 0.02092209133585608*I,
		-0.2877986626571256 - 0.04027358665471191*I,
		43.28298771155470174 + 0.01123031052443743*I}; //Compute these values in double precision with an external program (such as Mathematica).
	int i, nk = sizeof(ks)/sizeof(ks[0]);
	if (nk != sizeof(logdet_expc)/sizeof(logdet_expc[0])) {//Check the number of tests.
		fprintf(stderr, "[ERROR] Inconsistent number of tests, expected %d tests, aborting...\n", nk);
		return;
	}
	dcomplex* matrix = (dcomplex*)calloc(med->n*med->n, sizeof(dcomplex));
	dcomplex logdet;
	double err_real;
	for (i = 0; i < nk; i++) {//Loop on the test.
		build_ms_matrix(med, ks[i], matrix);
		logdet = tracelog_lu(med->n, matrix);
		err_real = fabs(creal(logdet)/creal(logdet_expc[i]) - 1.);
		printf("[INFO] Test #%2d | k=%5g%+5gi, log(det(M))/N=%10g%+10gi, expected=%10g%+10gi, err_real=%g \n",
			i, creal(ks[i]), cimag(ks[i]), creal(logdet), cimag(logdet), creal(logdet_expc[i]), cimag(logdet_expc[i]), err_real);
	}
	del_medium(med);
	free(med);
	free(matrix);
}

/**
 * Test the logarithmic derivative of the determinant of the multi-scattering matrix, (d/dk)log(det(M)) = Tr(M^-1 dM/dk).
 */
void test_dk_log_det_ms() {
	print_test("dk_log_det_ms()");
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {
		"dimension=2",
		"natom=50",
		"model=hardsphere 0.4",
		"shape=ball",
		"ratio=1.0"
	};
	uint64_t seed = 1; //Medium seed.
	
	parse_medium(med, sizeof(args)/sizeof(args[0]), args);
	print_param_medium(med);
	fill_medium(med, seed);
	//save_medium_points(med); exit(EXIT_SUCCESS);
	
	//Test several values of "k" against expected values of (d/dk)log(det(M)) computed in Wolfram Mathematica with the above settings for seed=1.
	dcomplex ks[] = {1. - 1.*I, 10. - 0.5*I, 10. - 10.*I, 0.1, 1., 3., 10., 100.};
	dcomplex res_expc[] = {36.41350253241055 + 75.48406981208467*I, -0.1605021020139895 + 182.8207081277391*I,
		13.33965080217422 + 247.5200299071803*I, 1.7999155479459534+8.217133651244348*I,
		14.404082456743884+25.986250874681836*I, 2.92720265720458+10.208149711288032*I,
		1.174691921299767+20.245454286532578*I, 341.9305967608437+19.998737802106373*I}; //Maintain this with an external program.
	int i, nk = sizeof(ks)/sizeof(ks[0]);
	if (nk != sizeof(res_expc)/sizeof(res_expc[0])) {//Check the number of tests.
		fprintf(stderr, "[ERROR] Inconsistent number of tests, expected %d tests, aborting...\n", nk);
		return;
	}
	dcomplex res;
	for (i = 0; i < nk; i++) {//Loop on the test.
		res = dk_log_det_ms(med, ks[i]);
		printf("[INFO] Test #%2d | k=%5g%+5gi, (d/dk)log(det(M))=%10g%+10gi, expected=%10g%+10gi, err=%g \n",
			i, creal(ks[i]), cimag(ks[i]), creal(res), cimag(res), creal(res_expc[i]), cimag(res_expc[i]), cabs(res/res_expc[i] - 1.));
	}
	del_medium(med);
	free(med);
}

/**
 * Test the inverse power iteration to find the eigenvector associated the smallest eigenvalue of a given matrix in absolute value.
 */
void test_inverse_iteration() {
	print_test("inverse_iteration()");
	dcomplex a[] = {-1.997 + 0.45*I, 1.875 - 1.47*I, -1.641 + 0.204*I, -1.74 - 1.613*I, -1.197 - 1.071*I, 
		1.883 - 1.765*I, 1.958 + 0.195*I, 1.278 + 1.105*I, -1.08 + 1.842*I};
	int n = (int)floor(sqrt(sizeof(a)/sizeof(a[0])));
	printf("[INFO] The matrix A is %dx%d.\n", n, n);
	
	dcomplex vec[n];
	constant_unit_cvector(n, vec); //Set "vec" to a normalized constant vector.
	dcomplex lambda;
	int maxit = 30;
	double toler = 1e-12;
	int verb = 1;  //Increase verbosity level for tests.
	
	inverse_iteration(n, a, vec, &lambda, maxit, toler, verb);
	dcomplex expc_vec[] = {0.16644079736120876, 0.4293668879243407 + 0.46904982189636174*I, 0.2043459505806509 + 0.7253802682650067*I};
	dcomplex expc_lambda = -0.38598154234503224 + 0.15816377074842095*I;
	
	print_cvector(n, vec, "Eigvec");
	print_cvector(n, expc_vec, "Expect");
	printf("[INFO] Error  = %g \n", cdistance(n, vec, expc_vec));
	
	printf("[INFO] Lambda = %.15g%+.15gi \n", creal(lambda), cimag(lambda));
	printf("[INFO] Expect = %.15g%+.15gi \n", creal(expc_lambda), cimag(expc_lambda));
	printf("[INFO] Error  = %g \n", cabs(lambda - expc_lambda));
}

/**
 * Test the computation of the trace of the inverse of a complex matrix, i.e., Tr(M^-1).
 */
void test_traceinv() {
	print_test("traceinv()");
	dcomplex a[] = {-1.997 + 0.45*I, 1.875 - 1.47*I, -1.641 + 0.204*I, -1.74 - 1.613*I, -1.197 - 1.071*I, 
		1.883 - 1.765*I, 1.958 + 0.195*I, 1.278 + 1.105*I, -1.08 + 1.842*I};
	int n = (int)floor(sqrt(sizeof(a)/sizeof(a[0])));
	printf("[INFO] The matrix A is %dx%d.\n", n, n);
	
	dcomplex trinv = traceinv(n, a);
	dcomplex trinv_expc = -2.5885257744954697 - 0.8796395991818305*I;
	
	printf("[INFO] Tr(A^-1) = %.15g%+.15gi \n", creal(trinv), cimag(trinv));
	printf("[INFO] Expected = %.15g%+.15gi \n", creal(trinv_expc), cimag(trinv_expc));
	printf("[INFO] Error    = %g \n", cabs(trinv - trinv_expc));
}

/**
 * Test the computation of the total cross section from the optical theorem: sigma(k) = -Im[phi·M^-1·phi]/k.
 */
void test_total_cross_section() {
	print_test("total_cross_section()");
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {
		"dimension=2",
		"natom=20",
		"model=hardsphere 0.4",
		"shape=ball",
		"ratio=1.0"
	};
	uint64_t seed = 1; //Medium seed used for the test.
	
	parse_medium(med, sizeof(args)/sizeof(args[0]), args);
	print_param_medium(med);
	fill_medium(med, seed);
	//save_medium_points(med); exit(EXIT_SUCCESS);
	
	double ks[] = {0.1, 1., 3., 10., 100.};
	double cs_expc[] = {20.52583442320853, 11.953344224711243, 9.862454676777112, 4.816530470222862, 0.0028798604611001526}; //Compute the expected total cross sections in double precision with an external program.
	
	int i, nk = sizeof(ks)/sizeof(ks[0]);
	if (nk != sizeof(cs_expc)/sizeof(cs_expc[0])) {//Check the number of tests.
		fprintf(stderr, "[ERROR] Inconsistent number of tests, expected %d tests, aborting...\n", nk);
		return;
	}
	double cs, err;
	for (i = 0; i < nk; i++) {//Loop on the tests.
		cs = total_cross_section(med, ks[i]);
		err = fabs(cs/cs_expc[i] - 1.); //Relative error.
		printf("[INFO] Test #%2d | k=%5g, sigma=%10g, expected=%10g, err=%g \n", i, ks[i], cs, cs_expc[i], err);
	}
	del_medium(med);
	free(med);
}

/**
 * Test the computation of the total scattering amplitude T(k,theta) in a specific direction "theta" in the (1x,1y) plane, assuming the particle is incident in the 1x direction.
 * This function also tests the generation of the plane wave vector.
 */
void test_scattering_amplitude() {
	print_test("Scattering amplitude");
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {
		"dimension=2",
		"natom=20",
		"model=hardsphere 0.4",
		"shape=ball",
		"ratio=1.0"
	};
	uint64_t seed = 1; //Medium seed used for the test.
	
	parse_medium(med, sizeof(args)/sizeof(args[0]), args);
	print_param_medium(med);
	fill_medium(med, seed);
	//save_medium_points(med); exit(EXIT_SUCCESS);
	
	double ks = 1.; //Fixed value of the wave number used for this test.
	double prefac = pow(ks/(2*PI), med->d-3)/(16*PI*PI); //Prefactor used to computes the differential cross section from the scattering amplitudes.
	double thetas[] = {0., 30, 90, 120, 180}; //Outgoing direction angle used for the tests (in degrees, interval 0:180).
	dcomplex atot_expc[] = {3.9719373901881365-11.9533442247113*I, 0.02483484083147164-8.152614522250493*I,
		-3.569044782687082+3.1417554077424166*I, -1.6098304636914875+4.718906749770716*I, 4.510670835267466+5.531912234005656*I}; //Compute the total scattering amplitude in double precision with an external program.
	double dcs_expc[] = {6.312830078587668, 2.6445877794572032, 0.8995718969592692, 0.9891334498146384, 2.027164640944822}; //Compute the differential cross section in double precision with an external program.
	
	int i, nth = sizeof(thetas)/sizeof(thetas[0]);
	if (nth != sizeof(atot_expc)/sizeof(atot_expc[0]) || nth != sizeof(dcs_expc)/sizeof(dcs_expc[0])) {//Check the number of tests.
		fprintf(stderr, "[ERROR] Inconsistent number of tests, expected %d tests, aborting...\n", nth);
		return;
	}
	int n = med->n, nrhs = 1, info; //Prepare for LAPACK's computation.
	dcomplex* matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));  //Allocate the MS matrix (OMP-Private).
	int* ipiv = (int*)calloc(n, sizeof(int));  //Pivot index array used by LAPACK (OMP-Private).
	dcomplex* a = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of scattering amplitudes (OMP-Private).
	dcomplex* b = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of output plane wave vector (OMP-Private).
	
	build_ms_matrix(med, ks, matrix); //Builds the matrix with the given value of "k".
	build_plane_wave(med, ks, 0., a); //Builds the plane wave oriented in the 1x direction (theta=0).
	zgesv_(&n, &nrhs, matrix, &n, ipiv, a, &n, &info); //Solve the linear system M*a = phi, so that "a" are now the scattering amplitudes per atom.
	
	dcomplex atot;
	double dcs, err_atot, err_dcs;
	for (i = 0; i < nth; i++) {//Loop on the scattering angles, one per bin.
		build_plane_wave(med, ks, thetas[i], b); //Builds the plane wave oriented in the output direction "th".
		atot = cdot_product(n, b, a); //Deduce the total scattering amplitude.
		dcs = prefac*(creal(atot)*creal(atot) + cimag(atot)*cimag(atot)); //Computes the differential cross section.
		err_atot = cabs(atot/atot_expc[i] - 1.); //Relative error on the total scattering amplitude.
		err_dcs  = fabs(dcs/dcs_expc[i] - 1.); //Relative error on the differential cross section.
		printf("[INFO] Test #%2d | th=%10g, atot=%10g%+10gi, expected=%10g%+10gi (err=%10g) | dcs=%10g, expected=%10g (err=%10g) \n",
			i, thetas[i], creal(atot), cimag(atot), creal(atot_expc[i]), cimag(atot_expc[i]), err_atot, dcs, dcs_expc[i], err_dcs);
	}
	del_medium(med);
	free(med);
	free(matrix);
	free(ipiv);
	free(a); free(b);
}

/**
 * Main function of the tests.
 */
int main(int argc, char** argv) {
	
	test_tracelog();
	test_dk_log_det_ms();
	test_inverse_iteration();
	test_traceinv();
	test_total_cross_section();
	test_scattering_amplitude();
	
	return 0;
}
