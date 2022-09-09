/****
 * @date Created on 2021-04-18 at 16:12:28 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing root-finding methods for complex functions, and fitting algorithms for real functions.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "real_vector_util.h"
#include "complex_vector_util.h"
#include "domain_util.h"

#define DELTA    1.e-4    //Step value used in finite difference. Too small steps lead to more cancellation errors.

/**********************************
 * BASIC LINEAR ALGEBRA FUNCTIONS
 *********************************/
/**
 * LAPACK's routine to solve over-determined linear systems of equations:
 */
void dgels_(char* trans, int* m, int* n, int* nrhs, double* a, int* lda, double* b, int* ldb, double* work, int* lwork, int* info);
void zgels_(char* trans, int* m, int* n, int* nrhs, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* work, int* lwork, int* info);

/**
 * Solves the given overdetermined real system of equations, A*X = B, where the real matrix "A" has "m" rows and "n" columns, and "B" has "m" rows.
 * In principle, we must have "m" > "n", and ideally "m" >> "n" for good fittings. This function interfaces with LAPACK's routine dgels_() for commodity.
 * @param "m" Number of rows of the matrices "a" and "b".
 * @param "n" Number of columns of the matrix "a".
 * @param "a" Least square rectangular real matrix, which is stored as an array in the LAPACK column-major format.
 * @param "b" Independent vector of the least square system (size "m"). On exit, the first "n" entries become the sought least square solution.
 * @returns "info" LAPACK's returned flag.
 */
int solve_least_square(int m, int n, double* a, double* b) {
	int info, nrhs = 1, lda = m, ldb = m, lwork = 2*n;
	double* work = (double*)calloc(lwork, sizeof(double)); //Allocate space for LAPACK's internal.
	dgels_("Notrans", &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
	if (info > 0) {
		fprintf(stderr, "[LAPACK] Element (%d,%d) of the triangular factor of A is zero. The least squares solution could not be computed.\n", info, info); fflush(stderr);
	}
	free(work);
	return info;
}

/**
 * Solves the given overdetermined complex system of equations, A*X = B, where the complex matrix "A" has "m" rows and "n" columns, and "B" has "m" rows.
 * In principle, we must have "m" > "n", and ideally "m" >> "n" for good fittings. This function interfaces with LAPACK's routine zgels_() for convenience.
 * @param "m" Number of rows of the matrices "a" and "b".
 * @param "n" Number of columns of the matrix "a".
 * @param "a" Least square rectangular complex matrix, which is stored as an array in the LAPACK column-major format.
 * @param "b" Independent vector of the least square system (size "m"). On exit, the first "n" entries become the sought least square solution.
 * @returns "info" LAPACK's returned flag.
 */
int solve_cleast_square(int m, int n, dcomplex* a, dcomplex* b) {
	int info, nrhs = 1, lda = m, ldb = m, lwork = 2*n;
	dcomplex* work = (dcomplex*)calloc(lwork, sizeof(dcomplex)); //Allocate space for LAPACK's internal.
	zgels_("Notrans", &m, &n, &nrhs, a, &lda, b, &ldb, work, &lwork, &info);
	if (info > 0) {
		fprintf(stderr, "[LAPACK] Element (%d,%d) of the triangular factor of A is zero. The least squares solution could not be computed.\n", info, info); fflush(stderr);
	}
	free(work);
	return info;
}

/**
 * Evaluates the given complex polynomial "poly" at the point "z" using backward evaluation (see Numerical Recipes, 3rd ed., chap. 5, p.202).
 * @param "d" Degree of the fitted polynomial, such that the highest power is x^d. The length of the list "poly" is then "d+1".
 * @param "poly" On exit, list of the coefficients of the polynomial in the standard order [a_0, a_1, a_2, ..., a_(d-1), a_d].
 * @param "z" Complex argument at which the polynomial is computed.
 * @return The value of the polynomial at "z".
 */
dcomplex eval_cpoly(int d, dcomplex* poly, dcomplex z) {
	if (d < 0) return 0.; //If the degree is negative, then returns zero.
	dcomplex res = poly[d];
	int i;
	for (i = d-1; i >= 0; i--) {//Backward evaluation is faster and more stable.
		res = res*z + poly[i];
	}
	return res;
}

/***********************************
 * SIMPLE LINEAR FITTING FUNCTIONS
 **********************************/
/**
 * Fit the real polynomial "poly" to the data points (xs[i], fs[i]).
 * @param "ns" Number of sample points in the "xs" and "fs" arrays.
 * @param "xs" List of the arguments of the sample points (size "ns").
 * @param "fs" List of the values of the sample points (size "ns").
 * @param "d" Degree of the fitted polynomial, such that the highest power is x^d. The length of the list "poly" is then "d+1".
 * @param "poly" On exit, list of the coefficients of the polynomial in the standard order [a_0, a_1, a_2, ..., a_(d-1), a_d].
 */
void fit_poly(int ns, double* xs, double* fs, int d, double* poly) {
	int i, j, np = d+1; //Number of polynomial coefficients, size of the polynomial array.
	if (ns < np) {
		fprintf(stderr, "[ERROR] Cannot fit a degree %d polynomial with only %d samples, needed at least %d, aborting...\n", d, ns, np); fflush(stderr);
		return;
	}
	double* a = (double*)calloc(ns*np, sizeof(double)); //Allocate space for the system matrix (column-major ordering).
	double* b = (double*)calloc(ns, sizeof(double)); //Allocate space for the independent term.
	double xp;
	for (i = 0; i < ns; i++) {//Loop on the rows of "a", the samples (xs, fs).
		xp = 1.;
		for (j = 0; j <= d; j++) {//Loop on the columns of "b", the polynomial coefficients (fitting parameters).
			a[i + j*ns] = xp; //Compute x[i]^j recursively.
			xp *= xs[i]; 
		}
		b[i] = fs[i];
	}
	solve_least_square(ns, np, a, b);
	copy_vector(np, b, poly);
	free(a);
	free(b);
}

/**
 * Fit the complex polynomial "poly" to the data points (zs[i], fs[i]).
 * @param "ns" Number of sample points in the "zs" and "fs" arrays.
 * @param "zs" List of the arguments of the sample points (size "ns").
 * @param "fs" List of the values of the sample points (size "ns").
 * @param "d" Degree of the fitted polynomial, such that the highest power is z^d. The length of the list "poly" is then "d+1".
 * @param "poly" On exit, list of the coefficients of the polynomial in the standard order [a_0, a_1, a_2, ..., a_(d-1), a_d].
 */
void fit_cpoly(int ns, dcomplex* zs, dcomplex* fs, int d, dcomplex* poly) {
	int i, j, np = d+1; //Number of polynomial coefficients, size of the polynomial array.
	if (ns < np) {
		fprintf(stderr, "[ERROR] Cannot fit a degree %d polynomial with only %d samples, needed at least %d, aborting...\n", d, ns, np); fflush(stderr);
		return;
	}
	dcomplex* a = (dcomplex*)calloc(ns*np, sizeof(dcomplex)); //Allocate space for the system matrix (column-major ordering).
	dcomplex* b = (dcomplex*)calloc(ns, sizeof(dcomplex)); //Allocate space for the independent term.
	dcomplex zp;
	for (i = 0; i < ns; i++) {//Loop on the rows of "a", the samples (zs, fs).
		zp = 1.;
		for (j = 0; j <= d; j++) {//Loop on the columns of "b", the polynomial coefficients (fitting parameters).
			a[i + j*ns] = zp; //Compute z[i]^j recursively.
			zp *= zs[i];
		}
		b[i] = fs[i];
	}
	solve_cleast_square(ns, np, a, b);
	copy_cvector(np, b, poly);
	free(a);
	free(b);
}

/**
 * Fit the complex polynomial "poly" to "ns" points of the given function "f(z)" computed in the domain "dom".
 * The points are picked up in the domain "dom" from the quasi-random (2,3) Halton sequence.
 * @param "f(z)" Complex function to be approximated by the polynomial "poly".
 * @param "dom" Domain of the complex plane in which the sample points (z, f(z)) are computed.
 * @param "ns" Number of sample points of "f(z)" to be computed.
 * @param "d" Degree of the fitted polynomial, such that the highest power is z^d. The length of the list "poly" is then "d+1".
 * @param "poly" On exit, list of the coefficients of the polynomial in the standard order [a_0, a_1, a_2, ..., a_(d-1), a_d].
 */
void fit_cpoly_from_func(dcomplex (*f)(dcomplex), Domain* dom, int ns, int d, dcomplex* poly) {
	dcomplex* zs = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the arguments of the sample points in the complex domain.
	dcomplex* fs = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the values of the function at the sample points, in the same order.
	set_halton_points(dom, ns, zs); //Pick up quasi-random uniform samples of the function in the domain.
	int i;
	for (i = 0; i < ns; i++) {//Loop on the sample points.
		fs[i] = f(zs[i]);  //Compute the function at the sample points zs[i].
	}
	fit_cpoly(ns, zs, fs, d, poly); //Fit the samples with polynomial by solving the linear system (calling LAPACK).
	free(zs);
	free(fs);
}

/**
 * Fit the complex Pade approximant "pade" to the data points (zs[i], fs[i]).
 * @param "ns" Number of sample points in the "z" and "f" arrays. The number of samples must be larger or equal to the number of coefficients "2*d+1".
 * @param "zs" List of the arguments of the sample points (size "ns").
 * @param "fs" List of the values of the sample points (size "ns").
 * @param "d" Degree of the balanced Pade approximant is [d/d]. So, the highest power is z^d in the numerator and the denominator.
 * @param "pade" List of the coefficients of the Pade approximant in the order [a_0, a_1, a_2, ..., a_d; b_1, b_2, ..., b_d], with b_0=1 by convention (size "2*d+1").
 */
void fit_cpade(int ns, dcomplex* zs, dcomplex* fs, int d, dcomplex* pade) {
	int i, j, np = 2*d + 1;  //Total number of Pade coefficients to be fitted.
	if (ns < np) {
		fprintf(stderr, "[ERROR] Cannot fit a [%d/%d] Pade with only %d samples, needed at least %d, aborting...\n", d, d, ns, np); fflush(stderr);
		return;
	}
	dcomplex* a = (dcomplex*)calloc(ns*np, sizeof(dcomplex));  //Allocate space for the system matrix (column-major ordering).
	dcomplex* b = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the independent term.
	dcomplex elem;  //Current matrix element.
	for (i = 0; i < ns; i++) {//Loop on rows of the matrix (samples).
		elem = 1.;  //Matrix element.
		for (j = 0; j < np; j++) {//Loop on columns of the matrix (parameters).
			a[i + j*ns] = elem;
			if (j == d) //Coefficients of the denominator.
				elem = -fs[i];
			elem *= zs[i];
		}
		b[i] = fs[i]; //Initializes "pade" as the independent vector term.
	}
	solve_cleast_square(ns, np, a, b);
	copy_cvector(np, b, pade); //Copy the result of "b" to the Pade array.
	free(a);
	free(b);
}

/**
 * Fit the complex Pade approximant "pade" of order [d/d] to "ns" points of the given function "f(z)" computed in the domain "dom".
 * The points are picked up in the domain "dom" from the quasi-random (2,3) Halton sequence.
 * @param "f(z)" Complex function to be approximated by "pade".
 * @param "dom" Domain of the complex plane in which the sample points (z, f(z)) are computed.
 * @param "ns" Number of sample points of "f(z)" to be computed.
 * @param "d" Degree of the balanced Pade approximant is [d/d]. So, the highest power is z^d in the numerator and the denominator.
 * @param "pade" List of the coefficients of the Pade approximant in the order [a_0, a_1, a_2, ..., a_d; b_1, b_2, ..., b_d], with b_0=1 by convention (size "2*d+1").
 */
void fit_cpade_from_func(dcomplex (*f)(dcomplex), Domain* dom, int ns, int d, dcomplex* pade) {
	dcomplex* zs = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the arguments of the sample points in the complex domain.
	dcomplex* fs = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the values of the function at the sample points, in the same order.
	set_halton_points(dom, ns, zs); //Pick up quasi-random uniform samples of the function in the domain.
	int i;
	for (i = 0; i < ns; i++) {//Loop on the sample points.
		fs[i] = f(zs[i]);  //Compute the function at the sample points zs[i].
	}
	fit_cpade(ns, zs, fs, d, pade); //Fit the samples with Pade approximants by solving the linear system (calling LAPACK).
	free(zs);
	free(fs);
}

/**************************
 * ROOT-FINDING FUNCTIONS 
 *************************/
/**
 * Declare LAPACK routine to find the eigenvalues of complex non-symmetric matrices.
 */
void zgeev_(char* jobvl, char* jobvr, int* n, dcomplex* a, int* lda, dcomplex* w, dcomplex* vl, int* ldvl, dcomplex* vr, int* ldvr, dcomplex* work, int* lwork, double* rwork, int* info);
/**
 * Finds all the roots of the given complex polynomial using the companion matrix method.
 * Note that this method is not guaranteed to provide machine accuracy for all the roots, especially for high-order polynomials.
 * In fact, finding the roots of a polynomial of larger degrees than about 100 makes little sense in double precision because of extreme loss of significance.
 * @param "d" Degree of the polynomial such that z^d is the highest power. So, the length of the array "poly" is just "d+1".
 * @param "poly" List of coefficients of the polynomial to be solved in standard order [a_0, a_1, a_2, ..., a_d] (size "d+1").
 * @param "roots" On exit, list of the computed roots of the given polynomials (size "d").
 */
void find_root_poly(int d, dcomplex* poly, dcomplex* roots) {
	dcomplex* a = (dcomplex*)calloc(d*d, sizeof(dcomplex));  //Allocate space for the matrix initially filled with zeros (column-major ordering).
	int i;
	for (i = 1; i < d; i++) {//Build the companion matrix.
		a[i + (d - 1)*d] = -poly[i]/poly[d];
		a[i + (i - 1)*d] = 1.;
	}
	a[(d - 1)*d] = -poly[0]/poly[d]; //Final matrix element.
	int lda = d, lwork = 4*d, info; //lwork=size of buffer array, can be arbitrary larger if necessary.
	dcomplex* work = (dcomplex*)calloc(lwork, sizeof(dcomplex)); //LAPACK Buffer.
	double* rwork = (double*)calloc(2*d, sizeof(double)); //LAPACK Buffer of size 2*d.
	zgeev_("Novec", "Novec", &d, a, &lda, roots, NULL, &d, NULL, &d, work, &lwork, rwork, &info); //Calls the zgeev_() subroutine to calculate eigenvalues.
	if (info > 0) {//Check for LAPACK error.
		fprintf(stderr, "[LAPACK] The QR algorithm failed to compute all the eigenvalues. Eigenvalues 1:%d have not converged (%.2f%%).\n", info, (100.0*info)/d); fflush(stderr);
	}
	sort_cvector(d, roots, creal); //Sort the roots in increasing order of real part.
	free(rwork);
	free(work);
	free(a);
}

/**
 * Finds a root of the given complex function "f(z)" using the damped secant method.
 * The damping factor "mu" is reduced as long as the new |f(z)| is larger than the current |f(z)|.
 * @param "f(z)" Double-precision complex function for which the roots are searched for.
 * @param "dom" Complex domain of root-finding. The iteration stops when going out of this domain.
 * @param "z" Initial guess of the root. On success, the argument "z" becomes a root of "f(z)".
 * @param "maxit" Maximum number of iterations (typically 30).
 * @param "toler" Tolerance over the relative error in the resulting "z" (typically 1e-12).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 * @returns 1 if the method has converged, 0 otherwise.
 */
int find_root_secant(dcomplex (*f)(dcomplex), Domain* dom, dcomplex* z, int maxit, double toler, int verb) {
	int s, conv = 0;  //Convergence flag which is 0 by default (no convergence).
	if (maxit <= 1 || toler <= 0.) {//Check for parameters.
		fprintf(stderr, "[ERROR] Invalid parameters (maxit=%d, toler=%g), aborting...\n", maxit, toler); fflush(stderr);
		return conv;
	}
	dcomplex za = (1. - DELTA*I)*(*z), zb = (1. + DELTA*I)*(*z);
	dcomplex fa = f(za), fb = f(zb), fn, zn, dz;
	double absfb = cabs(fb), absfn, absdz, tol;
	double mu, pow2;  //Damping factor.
	for (s = 0; s < maxit; s++) {//Secant iteration loop.
		if (fa == fb) {//Protects against division by zero.
			if (verb >= 1) {
				printf("[WARN] #%03d | Division by zero, aborting...\n", (s+1)); fflush(stdout);
			}
			break;
		}
		dz = fb*(zb - za)/(fb - fa); //Approximate f(z)/f'(z) by secant.
		absdz = cabs(dz);
		tol = toler*cabs(zb); //Absolute tolerance on abs(dz).
		if (verb >= 1) {
			printf("[STEP] #%03d | z=%20.14g%+20.14gi, abs(f)=%12g, abs(dz)=%12g \n", (s+1), creal(zb), cimag(zb), absfb, absdz); fflush(stdout);
		}
		if (!is_in_domain(dom, zb)) {
			if (verb >= 1) {
				printf("[WARN] #%03d | Root z=%20.14g%+20.14gi out of the domain, aborting...\n", (s+1), creal(zb), cimag(zb)); fflush(stdout);
			}
			break;
		}
		if (absdz < tol || absfb == 0.) {//Stopping criterion.
			*z = zb;  //Saves the root.
			conv = 1; //Sets the convergence flag to true.
			if (verb >= 1) {
				printf("[DONE] #%03d | Found root z=%20.14g%+20.14gi\n", (s+1), creal(zb), cimag(zb)); fflush(stdout);
			}
			break;
		}
		mu = 1.; pow2 = 2.;
		do { //The secant method is stabilized by damping (reduction to linear search).
			zn = zb - mu*dz;
			fn = f(zn); //Evaluates the time-costly function f(z) only here.
			absfn = cabs(fn);
			if (verb >= 2 && mu < 1.) {
				printf("[INFO] .... | z=%20.14g%+20.14gi, abs(f)=%12g, mu=%12g \n", creal(zn), cimag(zn), absfn, mu); fflush(stdout);
			}
			mu /= pow2; //Divides the damping factor by an increasingly large factor, mu = 2^(-n*(n+1)/2).
			pow2 *= 2.; //The reduction factor is increased in order to minimize the number of damping sub-steps.
		} while (absfn > absfb && mu*absdz > tol);
		za = zb; fa = fb; //Prepare for the next step.
		zb = zn; fb = fn; absfb = absfn;
	}
	return conv;
}

/**
 * Finds a root of the given complex function "f(z)" using the damped Muller method (local quadratic interpolation).
 * This method is more stable and requires less function evaluations that the secant method.
 * In particular, it is more suitable for the treatment of saddle points, f'(z)=0, than the secant method.
 * @param "f(z)" Double-precision complex function for which the roots are searched for.
 * @param "dom" Complex domain of root-finding. The iteration stops when going out of this domain.
 * @param "z" Initial guess of the root. On success, the argument "z" becomes a root of "f(z)".
 * @param "maxit" Maximum number of iterations (typically 30).
 * @param "toler" Tolerance over the relative error in the resulting "z" (typically 1e-12).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 * @returns 1 if the method has converged, 0 otherwise.
 */
int find_root_muller(dcomplex (*f)(dcomplex), Domain* dom, dcomplex* z, int maxit, double toler, int verb) {
	int s, conv = 0;  //Convergence flag which is 0 by default (no convergence).
	if (maxit <= 1 || toler <= 0.) {//Check for parameters.
		fprintf(stderr, "[ERROR] Invalid parameters (maxit=%d, toler=%g), aborting...\n", maxit, toler); fflush(stderr);
		return conv;
	}
	dcomplex za = (1. - DELTA*I)*(*z), zb = *z, zc = (1. + DELTA*I)*(*z);
	dcomplex fa = f(za), fb = f(zb), fc = f(zc), fn, zn, dz;
	dcomplex dfab, dfac, dfbc, d1f, d2f, delta, denom, denom1;
	double absfc = cabs(fc), absfn, absdz, tol;
	double mu, pow2;  //Damping factor.
	for (s = 0; s < maxit; s++) {//Muller iteration loop.
		dfab = (fa - fb)/(za - zb); //Computes the divided differences.
		dfac = (fa - fc)/(za - zc);
		dfbc = (fb - fc)/(zb - zc);
		d1f = dfbc + dfac - dfab; //First derivative, f'(zc). 
		d2f = (dfac - dfbc)/(za - zb); //Second derivative, f''(zc)/2.
		delta = csqrt(d1f*d1f - 4.*fc*d2f); //Solves the equation: fc + d1f*(z-zc) + d2f*(z-zc)^2 = 0, for "z".
		denom  = d1f + delta;
		denom1 = d1f - delta;
		if (cabs(denom1) > cabs(denom)) {//Chooses the largest denominator for stability.
			denom = denom1;
		}
		if (denom == 0. || denom != denom) {//Check for possible division by zeros, or NaN values.
			if (verb >= 1) {
				printf("[WARN] #%03d | Division by zero, aborting...\n", (s+1)); fflush(stdout);
			}
			break;
		}
		dz = 2.*fc/denom; //Step on "zc" such that zn = zc - dz.
		absdz = cabs(dz);
		tol = toler*cabs(zc); //Absolute tolerance on abs(dz).
		if (verb >= 1) {
			printf("[STEP] #%03d | z=%20.14g%+20.14gi, abs(f)=%12g, abs(dz)=%12g \n", (s+1), creal(zc), cimag(zc), absfc, absdz); fflush(stdout);
		}
		if (!is_in_domain(dom, zc)) {
			if (verb >= 1) {
				printf("[WARN] #%03d | Root z=%20.14g%+20.14gi out of the domain, aborting...\n", (s+1), creal(zc), cimag(zc)); fflush(stdout);
			}
			break;
		}
		if (absdz < tol || absfc == 0.) {//Stopping criterion.
			*z = zc;  //Saves the root.
			conv = 1; //Sets the convergence flag to true.
			if (verb >= 1) {
				printf("[DONE] #%03d | Found root z=%20.14g%+20.14gi\n", (s+1), creal(zc), cimag(zc)); fflush(stdout);
			}
			break;
		}
		mu = 1.; pow2 = 2.;
		do { //The Muller method is stabilized by damping, just like the secant method.
			zn = zc - mu*dz;
			fn = f(zn); //Evaluates the time-costly function f(z) only here.
			absfn = cabs(fn);
			if (verb >= 2 && mu < 1.) {
				printf("[INFO] .... | z=%20.14g%+20.14gi, abs(f)=%12g, mu=%12g \n", creal(zn), cimag(zn), absfn, mu); fflush(stdout);
			}
			mu /= pow2; //Divides the damping factor by an increasingly large factor, mu = 2^(-n*(n+1)/2).
			pow2 *= 2.; //The reduction factor is increased in order to minimize the number of damping sub-steps.
		} while (absfn > absfc && mu*absdz > tol);
		za = zb; fa = fb; //Prepare for the next step.
		zb = zc; fb = fc;
		zc = zn; fc = fn; absfc = absfn;
	}
	return conv;
}

/**
 * Approximates the locations of a certain number of roots of the given function "f(z)" within the given domain.
 * Use a balanced [n/n] Pade approximant fitted with points quasi-randomly picked up in the given domain.
 * @param "f(z)" Double-precision complex function.
 * @param "dom" Domain on which the roots are searched for.
 * @param "nz" Number of sought roots. On exit, this number becomes the actual number of approximate roots.
 * @param "z" Initial values are not used. On exit, computed estimates of the roots in the given domain "dom" (size pointed by "nz").
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 */
void guess_roots(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int verb) {
	int ns = 2*(*nz) + 1;     //Minimum required number of sample points to close the system of equations for the Pade approximant. The target number of zeros is "nz".
	dcomplex* pade = (dcomplex*)calloc(ns, sizeof(dcomplex));  //Allocate space for the Pade approximant, stored as [a0, a1, a2, ..., an, b1, b2, ..., bn], with b0=1 by convention.
	fit_cpade_from_func(f, dom, ns, *nz, pade); //Fit the samples with Pade approximants by solving the linear system (with LAPACK).
	find_root_poly(*nz, pade, z); //Root-find the numerator of the Pade approximant with the companion matrix method (with LAPACK).
	if (verb >= 1) {
		print_cvector(*nz, z, "Unfiltered roots"); fflush(stdout);
	}
	exclude_from_domain(dom, nz, z);  //Forgets about the roots outside the domain "dom" and reduce the logical size of "roots".
	free(pade);
}

/**
 * Finds several roots of the given complex function "f(z)" in the given domain "dom".
 * First guess the roots using the zeros of a fitted Pade approximant, then uses the Muller method separately on each root.
 * Note that the estimated number of zeros "nz" must be close to the actual number of zeros in the domain "dom".
 * Therefore, this method only makes sense if the number of roots "nz" is limited, and the function "f(z)" has not too many zeros in the domain "dom".
 * @param "f(z)" Double-precision complex function to be solved for "z".
 * @param "dom" Domain on which the roots are searched for.
 * @param "nz" Number of sought roots. On exit, this number becomes the number of successfully converged roots.
 * @param "z" On exit, computed roots which have successfully converged (size pointed by "nz").
 * @param "maxit" Maximum number of iterations of the Muller iterations (typically 30).
 * @param "toler" Tolerance on the relative error in the Muller iteration (typically 1e-12).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 */
void solve_from_guess(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int maxit, double toler, int verb) {
	guess_roots(f, dom, nz, z, verb); //First guess the roots using Pade approximant.
	int* flag = (int*)calloc(*nz, sizeof(int));
	int i;
	for (i = 0; i < *nz; i++) {//Loop on the roots to polish them with the Muller method.
		flag[i] = find_root_muller(f, dom, &z[i], maxit, toler, verb); //Stores whether the root has converged, or not.
	}
	filter_cvector(nz, z, flag);  //Filter non-converged roots.
	uniq_cvector( nz, z, toler);  //Remove duplicate roots within the prescribed tolerance on the relative error (note that this function does not assume pre-sorting).
	sort_cvector(*nz, z, creal);  //Sort the roots in ascending order of the real part.
	free(flag); //Free memory.
}

/**********************************************
 * MAEHLY-ABERTH-EHRLICH ROOT-FINDING METHOD
 *********************************************/
/**
 * Displays information about the current status of the Maehly-Aberth-Ehrlich root-finder.
 * Note that this function does not make any line return.
 * @param "s" Zero-based index of the current iteration.
 * @param "i" Zero-based index of the currently processed root.
 * @param "nz" Total number of complex roots.
 * @param "z" Current array of complex roots.
 * @param "flag" Current array of flags, one flag by root (0=in progress, 1=converged, -1=failed).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, >1=debug). This is also the number of displayed roots.
 */
void print_aberth_status(int s, int i, int nz, dcomplex* z, int* flag, int verb) {
	if (verb <= 0) return; //Printing something requires verbosity at least 1.
	int j, nprog = 0, nconv = 0, nfail = 0; //Counts the roots of the different category (in progress, converged, failed).
	for (j = 0; j < nz; j++) {//Loop on the roots.
		if (flag[j] == 0) {//Root still in progress.
			nprog++;
		}
		else if (flag[j] == 1) {//Converged root.
			nconv++;
		}
		else {//Failed root.
			nfail++;
		}
	}
	if (nprog != 0) {//If there are roots in progress, then show a summary.
		printf("\r[STEP] #%03d_%03d | nprog=%3d, nconv=%3d, nfail=%3d | sample: z[%3d]=%20.14g%+20.14gi ", (s+1), (i+1), nprog, nconv, nfail, i, creal(z[i]), cimag(z[i]));
	}
	else {
		printf("\r[DONE] #%03d_%03d | Found %3d/%3d roots,  nfail=%3d | sample: z[%3d]=%20.14g%+20.14gi ", (s+1), (i+1), nconv, nz, nfail, i, creal(z[i]), cimag(z[i]));
	}
	fflush(stdout);
}

/**
 * Returns 1 when the first "n" elements of the array "flag" are all non-zero.
 * Returns 0 if any of the elements of "flag" is equal to 0.
 */
int aberth_all_done(int n, int* flag) {//Private function.
	int i;
	for (i = 0; i < n; i++) {
		if (flag[i] == 0)
			return 0;
	}
	return 1;
}

/**
 * Finds at most "nz" distinct roots of the function "f(z)" using the Maehly-Aberth-Ehrlich iterative method
 * (aka Maehly's procedure in Numerical recipes) based on the given logarithmic derivative: dlogf(z) = f'(z)/f(z).
 * The roots "z" are not guessed by this function, but an external use of the Halton (2,3) sequence can be useful.
 * The method has the advantage of avoiding coincidences between the roots, so that every found root is unique.
 * Note that the Aberth method is guaranteed to converge as long as the function "f(z)" has no pole/singularity on the domain "dom" of interest.
 * To speedup the global convergence, this method can also suppress the exponential background of the form exp(a1*z + a2*z^2/2 + a3*z^3/3 + ...) from "f(z)".
 * This operation is implemented by fitting the polynomial a1 + a2*z + a3*z^2 + ... to "dlogf(z)" and then subtract it from "dlogf(z)" in the main iteration.
 * @param "dlogf(z)" Double-precision complex function equal to the logarithmic derivative f'(z)/f(z) of the function "f(z)" to be solved for "z".
 * @param "dom" Domain on which the roots are searched for.
 * @param "nz" Number of sought roots. On exit, this number becomes the number of successfully converged roots.
 * @param "z" Initial guess points of the roots. On exit, computed roots which have successfully converged (size pointed by "nz").
 * @param "nexpsup" Number of parameters involved in the suppression of the exponential background of f(z). Typical values are 2 or 4. Use nexpsup=0 to disable this feature.
 * @param "maxit" Maximum number of iterations (typically 30).
 * @param "toler" Tolerance on the relative error in the roots (typically 1e-12).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, >1=debug).
 */
void solve_aberth(dcomplex (*dlogf)(dcomplex), Domain* dom, int* nz, dcomplex* z, int nexpsup, int maxit, double toler, int verb) {
	uniq_cvector(nz, z, toler);      //Remove possibly duplicate starting points.
	exclude_from_domain(dom, nz, z); //Exclude starting points out of the domain.
	if (*nz <= 1) {//If there is less than 2 remaining guesses after filtering, then abort.
		fprintf(stderr, "[ERROR] Invalid initial roots (nz=%d), aborting...\n", *nz); fflush(stderr);
		return;
	}
	int s, i, j, nr = *nz; //Stores the original number of wanted roots.
	int* flag = (int*)calloc(nr, sizeof(int)); //Current status of each root (0=in progress, 1=converged, -1=failed).
	int d = nexpsup - 1;  //Degree of the polynomial used in the exponential suppression.
	dcomplex* poly;
	if (nexpsup > 0) {//If exponential suppression enabled, then fit a polynomial to dlogf(z), and subtract it in the main iterations.
		int ns = 6*nexpsup; //About 6 samples by fitting parameters should be enough.
		if (verb >= 1) {
			printf("[INFO] Exponential suppression enabled with nexpsup=%d. Pre-computing %d points...\n", nexpsup, ns); fflush(stdout);
		}
		poly = (dcomplex*)calloc(nexpsup, sizeof(dcomplex)); //Fitting polynomial approximating the background behavior of the function "dlogf(z)" to speedup convergence.
		fit_cpoly_from_func(dlogf, dom, ns, d, poly); //Computes extra samples from "dlogf(z)" to fit the polynomial.
	}
	else if (nexpsup == 0 && verb >= 1) {
		printf("[INFO] Exponential suppression disabled (nexpsup=%d).\n", nexpsup); fflush(stdout);
	}
	dcomplex invdz, dz, zn;
	for (s = 0; s < maxit; s++) {//Main loop of the Aberth root-finding method.
		for (i = 0; i < nr; i++) {//Loop on the roots.
			if (flag[i] == 0) {//Only process non-converged non-failed roots.
				invdz = dlogf(z[i]) - eval_cpoly(d, poly, z[i]); //Compute the time-costly function only here, and subtract the fitted polynomial to speedup convergence.
				for (j = 0; j < nr; j++) {//Loop on the roots to add the repulsion terms.
					if (j != i)
						invdz -= 1./(z[i] - z[j]);
				}
				dz = 1./invdz;
				zn = z[i] - dz;
				if (!is_in_domain(dom, zn)) {//If the new estimate is out of the domain, then the root estimate has failed.
					flag[i] = -1;
				}
				else if (cabs(dz) < toler*cabs(z[i])) {//Individual stopping criterion for each root. This may not detect exactly zero root, but nevertheless correct.
					flag[i] = 1;
				}
				else {//If nothing special happens, then accepts the new root estimate.
					z[i] = zn;
				}
				print_aberth_status(s, i, nr, z, flag, verb); //Prints the current status according to verbosity level.
			}
		}
		if (verb >= 1) {//Show a line break between iterations in verbose case.
			printf("\n");
		}
		if (aberth_all_done(nr, flag)) {//Stopping criterion: all roots must have converged or failed.
			break;
		}
	}
	filter_cvector(nz, z, flag); //Remove failed or non-converged roots by moving them to the end of the array "z", and reducing "nz".
	sort_cvector(*nz, z, creal); //Sort the roots in ascending order of the real part.
	free(flag);
}

/**
 * Finds at most "nz" distinct roots of the given function "f(z)" on the complex domain "dom" using the Maehly-Aberth-Ehrlich iterative method.
 * This function uses finite difference to evaluate the logarithmic derivative of the function "f(z)".
 * @deprecated Expected loss of speed for time-cotsly functions, but also some loss of significance,
 * possibly affecting the root-finding convergence speed. Used for tests only.
 */
void solve_aberth_num(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int maxit, double toler, int verb) {
	int nexpsup = 1; //Default value of the exponential suppression used for tests.
	dcomplex dlogf(dcomplex z) {//Defines the logarithmic derivative of f(z).
		return (f(z + DELTA) - f(z - DELTA))/(2*DELTA*f(z));
	};
	solve_aberth(dlogf, dom, nz, z, nexpsup, maxit, toler, verb);
}

/******************************
 * NON-LINEAR FITTING METHODS
 *****************************/
/**
 * Computes the residual, i.e., the sum of squared errors for the given model function "f(i,p)" with the given data "ys" and parameters "param".
 * In particular, this function computes sqrt((1/ns) Sum_i |y_i - f(i,p)|^2), where "ns" denotes the number of samples.
 */
double least_square_residual(double (*f)(int i, double* p), int ns, double* ys, double* param) {//Private function.
	double diff, r2 = 0.;
	int i;
	for (i = 0; i < ns; i++) {//Loop on the samples of data.
		diff = ys[i] - f(i, param);
		r2 += diff*diff;
	}
	return sqrt(r2/ns);
}

/**
 * Using the least-square system of equations for the given samples of data "ys" with the given parameters "param", find the correction on the parameters "dp".
 * Note that the Jacobian matrix is computed numerically using finite difference. Returns the info flag of LAPACK.
 */
int least_square_correction(double (*f)(int i, double* p), int ns, double* ys, int np, double* param, double* dp) {//Private function.
	double* a = (double*)calloc(np*ns, sizeof(double)); //Allocate space for the least-square "ns"-by-"np" Jacobian matrix.
	double* b = (double*)calloc(ns, sizeof(double)); //Allocate space for the independent term containing the "ns" residuals for all the samples.
	double* param_p = (double*)calloc(np*np, sizeof(double));  //Matrix containing the "np" parameter vectors modified by the finite difference in the "np" directions.
	double* param_m = (double*)calloc(np*np, sizeof(double));  //Matrix containing the "np" parameter vectors modified by the finite difference in the "np" directions.
	int i, j, info;
	for (i = 0; i < np; i++) {//Loop on the parameters, i.e., the rows of "param_p" and "param_m".
		for (j = 0; j < np; j++) {//Loop on the finite difference directions.
			param_p[i + j*np] = param[i];
			param_m[i + j*np] = param[i];
		}
		param_p[i + i*np] += DELTA; //Add or subtract the finite difference along the diagonal.
		param_m[i + i*np] -= DELTA;
	}
	for (i = 0; i < ns; i++) {//Loop on the samples, i.e., the rows of the Jacobian matrix.
		for (j = 0; j < np; j++) {//Loop on the parameters, i.e., the columns of the Jacobian matrix.
			a[i + j*ns] = (f(i,&param_p[j*np]) - f(i,&param_m[j*np]))/(2.*DELTA);
		}
		b[i] = ys[i] - f(i,param);
	}
	info = solve_least_square(ns, np, a, b); //Solve the least square system A*X = B using LAPACK.
	copy_vector(np, b, dp); //Copy the result into "dp".
	free(a); free(b);
	free(param_p);
	free(param_m);
	return info;
}
 
/**
 * Fit the parameters "param" of the given model function "f(i,param)" to the given samples of data "ys" using the Damped Gauss-Newton Algorithm.
 * The damping factor is divided by two until the residual decreases with respect to the previous one. The Jacobian matrix is computed numerically using finite differences.
 * Note that the tolerance on the result cannot be controlled in this function, because the final number of valid decimals depends on the curvature of the residual at the minimum.
 * @param "f(i,p)" Model function which is fitted to the data. "i" is the index of the sample, and "p" denote the parameters.
 * @param "ns" Number of samples of data.
 * @param "ys" Samples of data. Observed values of the model function "f(i,p)". Logical size "ns".
 * @param "np" Number of fitting parameters in the model function "f(i,p)".
 * @param "param" Initial values of the fitting parameters. Better guess provide better results. On success, the values become the found fitted parameters.
 * @param "maxit" Maximum number of Gauss-Newton iterations (typically 30).
 * @param "verb" Verbosity level (0=quiet, 1=verbose, 2=debug).
 * @returns 1 if the Gauss-Newton iteration has converged, 0 otherwise.
 */
int find_fit_gna(double (*f)(int i, double* p), int ns, double* ys, int np, double* param, int maxit, int verb) {
	int s, info, conv = 0;  //Convergence flag which is 0 by default (no convergence).
	if (maxit <= 1) {//Check for parameters.
		fprintf(stderr, "[ERROR] Invalid parameter maxit=%d, aborting...\n", maxit); fflush(stderr);
		return conv;
	}
	if (ns < np) {
		fprintf(stderr, "[ERROR] Cannot fit %d parameters with only %d points, aborting...\n", np, ns); fflush(stderr);
		return conv;
	}
	double toler = 1.e-14; //Tolerance on the relative variation of the residual between two consecutive steps.
	double* param_new = (double*)calloc(np, sizeof(double));   //New parameters.
	double* dp = (double*)calloc(np, sizeof(double));  //Correction on the parameters.
	double res = least_square_residual(f, ns, ys, param);  //Current residual.
	if (verb >= 1) {
		printf("[INIT] #%03d | res=%.16g \n", 0, res); fflush(stdout);
	}
	double res_new, mu, pow2;  //New residual, and damping factor "mu".
	for (s = 0; s < maxit; s++) {//Main loop of the Gauss-Newton method.
		info = least_square_correction(f, ns, ys, np, param, dp);  //Computes the correction "dp" on the parameters.
		if (info > 0) {//In case of LAPACK error while computing, then breaks the loop.
			break;
		}
		mu = 1.; pow2 = 2.; //Initializes the damping factor "mu".
		do { //The Gauss-Newton iteration is stabilized by damping.
			copy_vector(np, param, param_new); //Copies the content of "param" into "param_new".
			add_vector(np, param_new, dp, mu); //Equivalent to param_new += mu*dp. So, we have param_new = param + mu*dp.
			res_new = least_square_residual(f, ns, ys, param_new);  //Computes the new residual with the corrected parameters.
			if (verb >= 2 && mu < 1.) {
				printf("[INFO] .... | res=%.16g, mu=%12g \n", res_new, mu); fflush(stdout);
			}
			mu /= pow2; //Divides the damping factor by an increasingly large factor, mu = 2^(-n*(n+1)/2).
			pow2 *= 2.; //The reduction factor is increased in order to minimize the number of damping sub-steps.
		} while (res_new > res && mu > toler);
		copy_vector(np, param_new, param); //Saves the changes to the parameters.
		if (verb >= 1) {
			printf("[STEP] #%03d | res=%.16g, norm(dp)=%12g \n", (s+1), res_new, norm(np, dp)); fflush(stdout);
		}
		if (fabs(res_new - res) < toler*fabs(res)) {//Stopping criterion. Not (norm(np, dp) < toler*norm(np, param)), because the final precision of the parameters cannot be known in advance.
			conv = 1; //Sets the convergence flag to true.
			if (verb >= 1) {
				printf("[DONE] #%03d | Found fitted parameters.\n", (s+1)); fflush(stdout);
			}
			break;
		}
		res = res_new; //Saves the new residual in the ancient one.
	}
	free(param_new);
	free(dp);
	return conv;
}
