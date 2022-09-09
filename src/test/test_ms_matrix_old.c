/****
 * @date Created on 2020-08-09 at 19:31:37 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Test of the code providing the multiple scattering (MS) matrix creation utilities.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <time.h>           /* Standard Library for Time Measurement */
#include <omp.h>            /* Import the OpenMP Library for Parallelization */
#include "medium_util.h"    /* Import the Random Medium Management Utilities */
#include "ms_matrix_util.h" /* Import the Complex Library and the Multi-Scattering Matrix Utilities */

/**
 * Checks element by element that the two given n x n complex matrices (in column major ordering) are equal.
 * Prints the possible differences to standard output. Returns the number of found differences.
 */
int assert_matrix_equal(int n, dcomplex* mat1, dcomplex* mat2) {
	double tol = 1.0e-15; //Defines the tolerance on the relative error.
	double err; //Relative error between two terms.
	int i, j;
	int ndiff = 0; //Number of differences.
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			err = cabs(mat1[i + j*n] - mat2[i + j*n])/cabs(mat1[i + j*n]);
			if (err > tol) {
				printf("[WARN] Difference found on element (i=%d, j=%d) with relative error %f\n", i, j, err);
				ndiff++;
			}
		}
	}
	return ndiff;
}

/**
 * Compares the parallelized and the serial building of the multi-scattering matrix.
 * Asserts that they are equal, and estimate the speedup of the parallelized version over the serial one.
 */
void test_build_matrix(Medium* med, dcomplex k) {
	dcomplex* matser = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex)); //List of the MS matrix elements in column major order.
	dcomplex* matpar = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex));
	double timeser1 = omp_get_wtime();
	clock_t timeser2c = clock();
	build_matrix_serial(med, k, matser);
	timeser1 = omp_get_wtime() - timeser1;
	double timeser2 = (double)(clock() - timeser2c)/CLOCKS_PER_SEC;
	double timepar = omp_get_wtime();
	build_matrix_paral1(med, k, matpar);
	timepar = omp_get_wtime() - timepar;
	clock_t timechkc = clock();
	int ndiff = assert_matrix_equal(med->n, matser, matpar);
	double timechk = (double)(clock() - timechkc)/CLOCKS_PER_SEC;
	printf("[INFO] Time_serial_OMP = %fs, Time_serial = %fs, Time_parallel = %fs (speedup = %f) with %d difference(s) found. Check time = %fs\n", timeser1, timeser2, timepar, timeser1/timepar, ndiff, timechk);
	free(matser);
	free(matpar);
}

/**
 * Compare the performance of LAPACK eigenvalues computation, and the build of the MS matrix.
 */
void compare_build_vs_lapack(Medium* med, dcomplex k) {
	int i, ntest = 100;  //Increase the number of tests to get better statistics.
	dcomplex* matrix = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex)); //List of the MS matrix elements in column major order.
	
	double timebuildser = omp_get_wtime();
	for (i = 0; i < ntest; i++) {
		build_matrix_serial(med, k, matrix);
	}
	timebuildser = (omp_get_wtime() - timebuildser)/ntest;
	
	dcomplex* mu = (dcomplex*)malloc(med->n*sizeof(dcomplex)); //Eigenvalues of the MS matrix.
	int lwork = 3*med->n;
	dcomplex* work = (dcomplex*)malloc(lwork*sizeof(dcomplex)); //LAPACK Buffer.
	double* rwork = (double*)malloc(2*med->n*sizeof(double)); //LAPACK Buffer.
	
	double timelapack = omp_get_wtime();
	for (i = 0; i < ntest; i++) {
		eigvals(med->n, matrix, mu);
	}
	timelapack = (omp_get_wtime() - timelapack)/ntest;
	
	printf("[INFO] Time_build_serial = %fs, Time_LAPACK = %fs.\n", timebuildser, timelapack);
	
	free(mu);
	free(work);
	free(rwork);
	free(matrix);
}

int main(int argc, char** argv) {
	
	int d = 2;
	int n;
	uint64_t seed = 1;
	double ratio = 1.0;
	double alpha = 1.0;
	
	if (argc > 1) {
		int value;
		if (sscanf(argv[1], "%d", &value) == 1) {
			n = value;
		}
		else {
			printf("[ERROR] Incorrect number of atom '%s', expecting a positive integer.\n", argv[1]);
			return 1;
		}
	}
	else {
		printf("[ERROR] No number of atoms found, aborting. Usage: %s n, where n is the number of atoms.\n", argv[0]);
		return 1;
	}
	if (n < 1) {
		printf("[ERROR] The number of atoms '%d' should be a strictly positive integer, aborting.\n", n);
		return 1;
	}
	if (n > 5000) {
		printf("[ERROR] The number of atoms '%d' should be smaller than 5000 to avoid overloading the memory.\n", n);
		return 1;
	}
	
	Medium* med = new_medium(d, n, ratio, alpha);
	fill_medium(med, seed);
	print_param_medium(med, 1); //Display the values of the parameters of the medium to standard output for information.
	
	dcomplex k = 1.0 + 3.0*I;
	detect_overflow(med, cimag(k));
	
	//test_build_matrix(med, k);
	compare_build_vs_lapack(med, k);
	
	del_medium(med);
	return 0;
}
