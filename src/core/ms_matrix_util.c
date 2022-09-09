/****
 * @date Created on 2020-08-09 at 12:23 CEST
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Provides the tool to build and manipulate the multiple scattering (MS) matrix.
 ***/
#include <stdlib.h>               /* Standard Library for Memory Allocation */
#include <stdio.h>                /* Standard Library for Input and Output */
#include <math.h>                 /* Standard Library of Mathematical Function */
#include "complex_vector_util.h"  /* Import the Complex Vector Utilities */
#include "green_bessel.h"         /* Import the Complex Library and the Green and Bessel functions */
#include "medium_util.h"          /* Import the Random Medium Management Utilities */

/**
 * Define some constant macros:
 */
#define PI      3.1415926535897932385  //Value of pi.

/*********************************************
 * MULTI-SCATTERING MATRIX BUILDING FUNCTIONS
 ********************************************/
/**
 * Fills the given matrix with the appropriate elements of the multi-scattering matrix (column-major ordering) at the given complex wave number "k".
 * The MS matrix is defined as: M(k) = F(k)^-1 - G(k).
 * The matrix is assumed to be already allocated, i.e., N^2 double complex elements.
 */
void build_ms_matrix(Medium* med, dcomplex k, dcomplex* matrix) {
	dcomplex diag = invf(med, k);  //Call the inverse scattering amplitude of the atoms (from medium utils).
	dcomplex elem;  //Buffer value of matrix elements.
	int i, j;
	for (j = 0; j < med->n; j++) { //Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) { //Loop over the rows in the upper triangle (matrix is symmetric).
			elem = -free_green(med->d, k, atom_distance(med, i, j));  //Off-diagonal term.
			matrix[i + j*med->n] = elem; //Upper triangle.
			matrix[j + i*med->n] = elem; //Lower triangle.
		}
		matrix[(med->n + 1)*j] = diag; //Diagonal term.
	}
}

/**
 * Fills the given matrix with the elements of the normlized MS matrix (column-major ordering) at the given complex wave number "k".
 * The normalized MS matrix is defined as N(k) = i - G(k)/I(k,0), and is thus independent from the single-atom scattering model.
 */
void build_normalized_matrix(Medium* med, dcomplex k, dcomplex* matrix) {
	dcomplex izero = free_green_imag(med->d, k, 0.);  //Imaginary part of the free Green function, I(k,r)=-Im[G(k,r)], for r=0.
	dcomplex elem;  //Buffer value of matrix elements.
	int i, j;
	for (j = 0; j < med->n; j++) { //Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) { //Loop over the rows in the upper triangle (matrix is symmetric).
			elem = -free_green(med->d, k, atom_distance(med, i, j))/izero;  //Off-diagonal term.
			matrix[i + j*med->n] = elem; //Upper triangle.
			matrix[j + i*med->n] = elem; //Lower triangle.
		}
		matrix[(med->n + 1)*j] = I; //Diagonal term.
	}
}

/**
 * Fills the given matrix with the elements of the k-derivative multi-scattering matrix (column-major ordering) at the given complex wave number "k".
 * The matrix is assumed to be already allocated, i.e., N^2 double complex elements.
 */
void build_dk_ms_matrix(Medium* med, dcomplex k, dcomplex* matrix) {
	const dcomplex dk = 1.e-4;  //Delta used for numerical derivative (not too small to avoid loss of significance).
	dcomplex diag = (invf(med, k + dk) - invf(med, k - dk))/(2*dk); //This is only an approximation to 7-8 decimal places, but enough for the Newton-Raphson method.
	dcomplex elem;  //Buffer value of matrix elements.
	int i, j;
	for (j = 0; j < med->n; j++) { //Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) { //Loop over the rows in the upper triangle (matrix is symmetric).
			elem = -dk_free_green(med->d, k, atom_distance(med, i, j));  //Off-diagonal term.
			matrix[i + j*med->n] = elem; //Upper triangle.
			matrix[j + i*med->n] = elem; //Lower triangle.
		}
		matrix[(med->n + 1)*j] = diag; //Diagonal term.
	}
}

/**
 * Fills the given vector "vec" with the values of the unit-amplitude plane wave of wave number "k" and angle "theta" (in degrees)
 * with respect to the 1x direction in the (1x,1y) plane.
 * The vector is assumed to be already allocated, i.e., with N double complex elements.
 */
void build_plane_wave(Medium* med, dcomplex k, double theta, dcomplex* vec) {
	int i, d = med->d, n = med->n;
	theta *= PI/180; //Converts from degrees to radians. "theta" is now measured in radians in the interval 0:pi.
	double costheta = cos(theta), sintheta = sin(theta);
	if (d == 1) {//In the 1D case, the direction is either theta=0 or theta=pi.
		if (costheta > 0.)
			costheta = +1.;
		else
			costheta = -1.;
		for (i = 0; i < n; i++) {//Loop on the atoms.
			vec[i] = cexp(I*k*costheta*med->pos[i]);
		}
	}
	else {//In any dimension d >= 2, the angle "theta" is taken into account.
		for (i = 0; i < n; i++) {//Loop on the atoms.
			vec[i] = cexp(I*k*(costheta*med->pos[i*d] + sintheta*med->pos[i*d+1]));
		}
	}
}

/**
 * Displays the first s x s block, i.e., the scope, of the given n x n complex matrix (in column-major ordering). Does not check for bounds.
 */
void print_matrix(int n, dcomplex* matrix, int s) {
	if (s > n) {
		s = n;
	}
	printf("[INFO] M = \n");
	int i, j;
	for (i = 0; i < s; i++) {
		printf("[ROW%d]", i+1);
		for (j = 0; j < s; j++) {
			printf("\t%f%+fi", creal(matrix[i + j*n]), cimag(matrix[i + j*n]));
		}
		if (s < n) {
			printf("\t...");
		}
		printf("\n");
	}
	fflush(stdout);
}

/****************************
 * LAPACK MATRIX OPERATIONS
 ***************************/
/**
 * Declare subroutines from the LAPACK library:
 * Reminder: z=Double complex (complex*16, totalling 16 bytes), ge=General matrix, ev=Eigenvalues, trf=Triangular factorization (LU factorization), sv=Linear solve.
 */
void zgeev_(char*, char*, int*, dcomplex*, int*, dcomplex*, dcomplex*, int*, dcomplex*, int*, dcomplex*, int*, double*, int*);  //Eigenvalues of complex general matrix.
void zgeqrf_(int* m, int* n, dcomplex* matrix, int* lda, dcomplex* tau, dcomplex* work, int* lwork, int* info);  //QR decomposition of a general complex matrix (improved stability).
void zgetrf_(int* m, int* n, dcomplex* matrix, int* lda, int* ipiv, int* info);  //LU decomposition of a complex general matrix.
void zgetrs_(char* trans, int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv, dcomplex* b, int* ldb, int* info); //Solves a system from the LU factorization given by zgetrf_().
void zgesv_(int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv, dcomplex* b, int* ldb, int* info);  //Solves a system using LU decomposition method.
void zgetri_(int* n, dcomplex* a, int* lda, int* ipiv, dcomplex* work, int* lwork, int* info); //Inverse matrix.

/**
 * Computes the N eigenvalues "mu" of the given general N x N complex matrix using LAPACK's zgeev_().
 */
void eigvals(int n, dcomplex* matrix, dcomplex* mu) {
	char nochar = 'N';  //'N'=Does not compute left/right eigenvectors.
	dcomplex* vl = NULL; //Empty eigenvectors.
	dcomplex* vr = NULL;
	int ldv = 1; //Leading dimension of eigenvectors (1 if not needed).
	int lwork = 3*n; //Size of buffer array, can be arbitrary larger if necessary.
	dcomplex* work = (dcomplex*)malloc(lwork*sizeof(dcomplex)); //LAPACK Buffer.
	double* rwork = (double*)malloc(2*n*sizeof(double)); //LAPACK Buffer.
	int info;
	zgeev_(&nochar, &nochar, &n, matrix, &n, mu, vl, &ldv, vr, &ldv, work, &lwork, rwork, &info); //Calls the zgeev_() subroutine to calculate eigenvalues.
	if (info > 0) {//Check for LAPACK error.
		printf("[LAPACK] The QR algorithm failed to compute all the eigenvalues. Eigenvalues 1:%d have not converged (%.2f%%).\n", info, (100.0*info)/n);
	}
	free(rwork);
	free(work);
}

/**
 * Computes the trace Tr[G⁺G]/N from the given multiple scattering matrix "M" of size "n".
 */
double tracegtwo(int n, dcomplex* matrix) {
	double tr = 0.;
	dcomplex g;
	int i, j;
	for (j = 0; j < n; j++) {//Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) {//Loop over the rows in the upper triangle (matrix is symmetric).
			g = matrix[i + j*n];
			tr += creal(g)*creal(g) + cimag(g)*cimag(g);
		}
	}
	return 2*tr/n; //Factor of two to account for the lower triangle.
}

/**
 * Computes the normalized trace-log function, i.e., Tr[ln(M)]/N, of the given general complex matrix "M" of size "n",
 * using the determinant method and LU decomposition (partial pivoting) with LAPACK's zgetrf_().
 * This method is the fastest and the most stable to far from the real "k" axis. It is especially accurate for large real parts in the complex "k" plane.
 */
dcomplex tracelog_lu(int n, dcomplex* matrix) {
	int i, info, ipiv[n], sign = 1;
	zgetrf_(&n, &n, matrix, &n, ipiv, &info); //LU decomposition with partial pivoting (most stable LAPACK factorization).
	if (info > 0) {//Check for LAPACK error.
		printf("[LAPACK] Element U(%d,%d) is exactly zero, the determinant will be singular.\n", info, info);
	}
	dcomplex trlog = 0.0; //Initializes the sought value of the trace-log function.
	for (i = 0; i < n; i++) {//The sum of logarithm prevents possible overflow if the diagonal elements were multiplied.
		trlog += clog(matrix[i+i*n]);
		if (ipiv[i] != i+1) {//Computes the permutation signature of the partial pivoting (even/odd number of binary interchanges):
			sign = -sign;
		}
	}
	if (sign == -1) {//If negative sign, the determinant changes sign, so its logarithm gets +i*pi. NB: Opposite sign -i*pi would be correct, too.
		trlog += I*PI;
	}
	return trlog/n; //Final normalization.
}

/**
 * Computes the normalized trace-log function, i.e., Tr[ln(M)]/N, of the given general complex matrix "M" of size "n",
 * using the spectral decomposition method with LAPACK's zgeev_().
 * @deprecated LAPACK's iterative QR algorithm is much slower than LU decomposition method below. Also less stable than LU far from the real "k" axis.
 */
dcomplex tracelog_ev(int n, dcomplex* matrix) {
	dcomplex* mu = (dcomplex*)malloc(n*sizeof(dcomplex)); //Allocate space for eigenvalues.
	eigvals(n, matrix, mu);
	dcomplex trlog = 0.0;
	int i;
	for (i = 0; i < n; i++) {//Loop over the eigenvalues of the matrix.
		trlog += clog(mu[i]);
	}
	free(mu); //Free eigenvalues array "mu".
	return trlog/n; //Final normalization.
}

/**
 * Computes the normalized trace-log function, i.e., Tr[ln(M)]/N, of the given general complex matrix "M" of size "n",
 * using the determinant method and QR decomposition (Householder reflections method) with LAPACK's zgeqrf_().
 * Note that the imaginary part of the function is irrelevant, because the sign of the determinant is not computed from Q.
 * @deprecated This version is less stable than the LU method far from the real "k" axis, especially for large real parts of "k".
 */
dcomplex tracelog_qr(int n, dcomplex* matrix) {
	int i, info, lwork = n;
	dcomplex tau[n], work[lwork], trlog = 0.0;
	zgeqrf_(&n, &n, matrix, &n, tau, work, &lwork, &info);  //QR decomposition.
	for (i = 0; i < n; i++) {//The sum of logarithm prevents possible overflow if the diagonal elements were multiplied.
		trlog += clog(matrix[i+i*n]);
	}
	return trlog/n; //Final normalization.
}

/**
 * Computes the logarithmic derivative of the determinant of the multiple scattering matrix, det(M), at the given complex wave number "k".
 * This function uses the formula (d/dk)ln(det(M)) = Tr(M^-1 dM/dk) to perform the calculation, hence avoiding the numerical evaluation of the derivative of the discontinuous function ln(det(M)).
 * This function is at the core of the Newton-Raphson algorithm used to find the roots of the resonance equation: det(M)=0. This function is also used to compute the density of states (DOS).
 */
dcomplex dk_log_det_ms(Medium* med, dcomplex k) {
	int n = med->n, lda = n, ldb = n, nrhs = n, info; //Prepare for the LAPACK calls.
	dcomplex* a = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	dcomplex* b = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	int* ipiv = (int*)calloc(n, sizeof(int)); //Allocate space for the permutation index array used by LAPACK.
	build_ms_matrix(med, k, a);     //Computes the matrix M(k) in "a".
	build_dk_ms_matrix(med, k, b);  //Computes the matrix dM/dk in "b".
	zgetrf_(&n, &n, a, &lda, ipiv, &info);  //Computes the LU decomposition of the matrix "a".
	if (info > 0) {//Check for possible LAPACK error.
		printf("[LAPACK] Element U(%d,%d) is exactly zero, the inverse matrix will be singular.\n", info, info);
	}
	zgetrs_("Notrans", &n, &nrhs, a, &lda, ipiv, b, &ldb, &info); //Computes the matrix M^-1*dM/dk overwriting "b".
	dcomplex trace = 0.;
	int i;
	for (i = 0; i < n; i++) {//Loop to compute the trace of "b".
		trace += b[(n+1)*i];
	}
	free(a); //Frees the memory.
	free(b);
	free(ipiv);
	return trace;
}

/**
 * Finds the eigenvector of the given "n"-by-"n" matrix "a" (column-major ordering) corresponding to the smallest eigenvalue of "a" in absolute value.
 * The resulting complex eigenvector is stored in "vec", and the corresponding eigenvalue is stored in "lambda".
 * Note that this function uses the elementary inverse power iteration and is efficient only
 * if the smallest eigenvalue of "a" is close to zero compared to the other eigenvalues. On exit, the matrix "a" is altered by LAPACK.
 * @param "n" Size of the matrix "a" and of the vectors.
 * @param "a" Square complex matrix for which the eigenvector is sought. On exit, this matrix is altered by LAPACK.
 * @param "vec" Initial guess of the sought eigenvector. On success, eigenvector corresponding to the smallest eigenvalue of "a".
 * @param "lambda" On success, smallest eigenvalue of "a" associated with the eigenvector "vec".
 * @param "maxit" Number of inverse power iterations, typically 50.
 * @param "toler" Tolerance on the relative distance between two successive steps of the inverse power iteration, typically 1e-10.
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 * @return 1 if the method has converged within the prescribed "maxit" iterations, 0 otherwise.
 */
int inverse_iteration(int n, dcomplex* a, dcomplex* vec, dcomplex* lambda, int maxit, double toler, int verb) {
	int lda = n, ldb = n, nrhs = 1, info;
	int* ipiv = (int*)calloc(n, sizeof(int)); //Allocate space for the permutation index array used by LAPACK.
	zgetrf_(&n, &n, a, &lda, ipiv, &info);  //Computes the LU decomposition of the matrix "a" only once since O(n^3).
	if (info > 0) {//Check for possible LAPACK error.
		printf("[LAPACK] Element U(%d,%d) is exactly zero, the inverse power iteration will fail.\n", info, info); fflush(stdout);
	}
	dcomplex* vold = (dcomplex*)calloc(n, sizeof(dcomplex)); //Allocate space for the previous eigenvector at each step.
	dcomplex fac = cnormalize(n, vec); //Normalizes the guess vector "vec", just in case, and saves the normalization factor.
	double delta;  //Distance between two consecutive eigenvectors.
	int s, conv = 0; //Convergence flag.
	for (s = 0; s < maxit; s++) {//Inverse power method iterations.
		copy_cvector(n, vec, vold);  //Copy "vec" into "vold" as a back up.
		zgetrs_("Notrans", &n, &nrhs, a, &lda, ipiv, vec, &ldb, &info); //Computes vec = A^-1*vold, only O(n^2).
		fac = cnormalize(n, vec); //Normalizes the new vector "vn", and computes the factor such that: original_vec = fac*normalized_vec.
		delta = cdistance(n, vec, vold);
		if (verb >= 1) {
			printf("[STEP] #%03d | fac=%18.12g%+18.12gi, delta=%g \n", s, creal(fac), cimag(fac), delta); fflush(stdout);
		}
		if (delta < toler) {//Stopping criterion. Note that cnormalize() also locks the complex phase of "vn".
			conv = 1;
			break;
		}
	}
	*lambda = 1./fac; //Compute the eigenvalue from the amplification factor.
	if (verb >= 1) {
		if (conv)
			printf("[INFO] Found eigenvector with eigenvalue lambda=%g%+gi\n", creal(*lambda), cimag(*lambda));
		else
			printf("[WARN] Eigenvector not found (last lambda=%g%+gi)\n", creal(*lambda), cimag(*lambda));
		fflush(stdout);
	}
	free(vold);
	free(ipiv);
	return conv;
}

/**
 * Computes the trace of the inverse of the given n-by-n complex matrix "a" in column-major format.
 */
dcomplex traceinv(int n, dcomplex* a) {
	int m = n, lda = n, lwork = 4*n, i, info;
	int* ipiv = (int*)calloc(n, sizeof(int));
	dcomplex* work = (dcomplex*)calloc(lwork, sizeof(dcomplex));
	zgetrf_(&m, &n, a, &lda, ipiv, &info);  //LU decomposition of a complex general matrix.
	if (info > 0) {
		printf("[LAPACK] Element U(%d,%d) is exactly zero, the inverse will be singular.\n", info, info); fflush(stdout);
	}
	zgetri_(&n, a, &lda, ipiv, work, &lwork, &info); //Inverse of the matrix.
	dcomplex tr = 0.;
	for (i = 0; i < n; i++) {//Compute the trace of the resulting matrix.
		tr += a[i + i*n];
	}
	free(ipiv);
	free(work);
	return tr;
}

/**
 * Computes the total cross section of the given medium "med" at the given wave number "k" in units of sp^(d-1).
 * The initial wave is assumed to be an incident plane wave oriented in the main (1x) direction.
 * The computation method comes from the optical theorem: sigma(k) = -Im[phi·M^-1·phi]/k.
 * The total cross section cannot be continued to the complex "k" plane without doubling the computational cost.
 * This is why the function only takes real-valued "k".
 */
double total_cross_section(Medium* med, double k) {
	int n = med->n, lda = n, ldb = n, nrhs = 1, ipiv[n], info; //Prepare for the LAPACK calls.
	dcomplex* matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	dcomplex* b = (dcomplex*)calloc(n, sizeof(dcomplex));
	dcomplex* phi = (dcomplex*)calloc(n, sizeof(dcomplex));
	build_ms_matrix(med, k, matrix);  //Computes the matrix M(k).
	int i;
	dcomplex z;
	for (i = 0; i < n; i++) {//Computes the incident plane wave on each atom.
		z = cexp(I*k*med->pos[i*med->d]); //Only the "x" component, since the plane wave is oriented in the 1x direction.
		b[i] = z; //This vector will be changed to the scattering amplitudes after use of LAPACK's routine.
		phi[i] = z; //This vector remains unaltered.
	}
	zgesv_(&n, &nrhs, matrix, &lda, ipiv, b, &ldb, &info);  //Solves the system using LU decomposition method.
	if (info > 0) {//Check for LAPACK errors.
		printf("[LAPACK] Element U(%d,%d) is exactly zero, the solution could not be computed.\n", info, info); fflush(stdout);
	}
	double cs = -cimag(cdot_product(n, phi, b))/k; //Computes the resulting total cross section.
	free(matrix); //Frees the memory.
	free(b);
	free(phi);
	return cs;
}
