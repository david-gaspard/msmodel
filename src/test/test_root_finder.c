/****
 * @date Created on 2021-04-18 at 18:11:38 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to test the root finding methods.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "tag_strings.h"
#include "assertion.h"
#include "real_vector_util.h"
#include "complex_vector_util.h"
#include "domain_util.h"
#include "root_finder.h"

#define TOLER      1.e-13     //General tolerance used for most tests on relative errors of numerical computations.
#define TWOPI      6.2831853071795864769

/*********************
 * PRELIMINARY TESTS: 
 ********************/
/**
 * Test of parametric function and variable scoping.
 */
int test_parametric_function() {
	print_test("parametric function");
	int nf = 0; //Number of failed tests.
	int delta = 0;
	int increment(int value) {
		return value + delta;
	}
	int x = 6;
	int res = increment(x);
	int res_expc = x+delta;
	nf += (res != res_expc);
	printf("[INFO] delta=%d, inc(%d)=%d, expected=%d \n", delta, x, res, res_expc);
	
	delta = 1; //Here modifies the parameter of the function increment(int).
	res = increment(x);
	res_expc = x+delta;
	nf += (res != res_expc);
	printf("[INFO] delta=%d, inc(%d)=%d, expected=%d \n", delta, x, res, res_expc);
	
	if (nf == 0) {
		printf(STR_PASS "Tests passed.\n");
	}
	return nf;
}

/**
 * Test the evulation of complex polynomials.
 */
int test_eval_cpoly() {
	print_test("eval_cpoly()");
	dcomplex poly[] = {-0.435 - 4.337*I, 3.682 + 2.83*I, 2.043 + 4.074*I, 2.95 - 4.19*I, -4.595 - 0.134*I,
		4.578 + 3.248*I, -4.916 + 1.845*I, -2.487 - 4.362*I, -4.857 - 4.137*I, 2.439 + 1.584*I};
	dcomplex z[] = {0, 1, -1, I, -I, 2+5*I, -1-3*I, -2-7*I, 4-6*I};
	dcomplex p_expc[] = {-0.435 - 4.337*I, -1.598 - 3.579*I, -23.922 - 1.799*I, -23.228 - 4.291*I, 9.2 - 24.763*I,
		6.903302368e6 - 7.434526259e6*I, -68407.251 + 113273.605*I, -1.63661892832e8 + 1.04628931181e8*I, -9.7142808275e7 - 1.06715262285e8*I};
	int d = sizeof(poly)/sizeof(poly[0]) - 1;
	int i, n = sizeof(z)/sizeof(z[0]);
	int n_expc = sizeof(p_expc)/sizeof(p_expc[0]);
	dcomplex p[n];
	for (i = 0; i < n; i++) {//Loop on the tests.
		p[i] = eval_cpoly(d, poly, z[i]);
	}
	print_cvector(n, p, "Eval results");
	print_cvector(n_expc, p_expc, "Expc results");
	return assert_cvector_close(n, p, n_expc, p_expc, TOLER);
}

/***********************************
 * TEST THE ROOT-FINDING UTILITIES:
 **********************************/

/**
 * Tests the real polynomial fitting function.
 */
int test_fit_poly() {
	print_test("fit_poly()");
	double xs[] = {0.635, 0.579, -0.517, 0.084, -0.208, -0.576, -0.154, 0.954, 0.851, -0.414};
	double fs[] = {2.442, 2.309,  2.746, 1.669,  2.529,  3.629,  1.969, 4.970, 4.016,  2.565};
	int ns = sizeof(xs)/sizeof(xs[0]);
	int d = 3, np = d+1;
	double poly[np];
	fit_poly(ns, xs, fs, d, poly);
	double poly_expc[] = {1.7580096943876664, -1.6099797279792991, 3.052116483868505, 2.2537615124971877};
	int n_expc = sizeof(poly_expc)/sizeof(poly_expc[0]);
	print_vector(np, poly, "Poly");
	print_vector(n_expc, poly_expc, "Expc");
	return assert_vector_close(np, poly, n_expc, poly_expc, TOLER);
}

/**
 * Test the complex Pade fitting function.
 */
int test_fit_cpade() {
	print_test("fit_cpade()");
	dcomplex zs[] = {0.635 - 0.777*I, -0.517 - 0.869*I, -0.208 + 0.401*I, -0.154 - 0.505*I, 
		0.851 + 0.156*I, 0.161 - 0.742*I, -0.219 + 0.64*I, 0.038 - 0.662*I, -0.976 - 0.366*I, -0.217 - 0.082*I};
	dcomplex fs[] = {0.231 - 0.802*I, 0.825 - 1.189*I, -0.014 + 2.407*I, 0.659 - 1.347*I, -0.486 - 0.182*I,
		-0.193 - 0.690*I, 0.259 + 1.457*I, -0.011 - 0.804*I, 1.751 - 0.972*I, -3.237 - 8.956*I};
	int ns = sizeof(zs)/sizeof(zs[0]);
	int d = 1; //Degree of the Pade [d/d].
	int np = 2*d + 1;  //Number of Pade parameters.
	dcomplex pade[np];
	fit_cpade(ns, zs, fs, d, pade);  //Fit the Pade approximant to the data.
	dcomplex pade_expc[] = {-2.6405533329787456 - 0.5025177225177638*I, 1.870390152106856 - 0.23376120382041643*I, 3.385043965133553 + 0.18079515867305634*I}; //Expected result.
	int n_expc = sizeof(pade_expc)/sizeof(pade_expc[0]);
	print_cvector(np, pade, "Pade");
	print_cvector(n_expc, pade_expc, "Expc");
	return assert_cvector_close(np, pade, n_expc, pade_expc, TOLER);
}

/**
 * Tests the function to find the roots of a polynomial.
 */
int test_find_root_poly() {
	print_test("find_root_poly()");
	dcomplex poly[] = {1.904 - 1.729*I, -2.331 + 1.492*I, 1.737 - 0.463*I, -1.873 - 1.515*I, -1.552 + 2.863*I,
		-2.606 + 1.951*I, 0.253 + 2.552*I, -1.613 + 0.468*I, -0.624 - 1.243*I, 1.203 - 1.752*I};
	int np = sizeof(poly)/sizeof(poly[0]);
	int d = np-1; //The polynomial degree is also the number of roots.
	dcomplex roots[d];
	find_root_poly(d, poly, roots);
	dcomplex roots_expc[] = {-1.3112405355944026 - 0.28719921917538255*I, -0.8594174347433288 + 0.8486268297218884*I, -0.838880801550193 - 0.5405345084565891*I, 
		-0.06477069346953607 - 0.9731203147377844*I, 0.0012397260242508897 + 1.3422148459497214*I, 0.25585646599966966 + 0.6948501698668312*I, 
		0.42218965542036163 - 0.678396304677164*I, 0.6825979571449342 - 0.04262208329559775*I, 1.3964737727470218 + 0.20929138462686797*I};
	print_cvector(d, roots, "Roots");
	print_cvector(d, roots_expc, "Expec");
	return assert_cvector_close(d, roots, sizeof(roots_expc)/sizeof(roots_expc[0]), roots_expc, TOLER);
}

/**
 * Tests the secant root-finder for the given initial values "z" (array of size "n"), and checks that all the resulting roots belong to the "z_expc" array.
 * Returns 1 on failure, otherwise zero.
 */
int test_find_root_method_from_input(char* solvername, int (*find_root)(dcomplex (*f)(dcomplex), Domain*, dcomplex*, int, double, int), char* funcname, dcomplex (*f)(dcomplex), Domain* dom, int n, dcomplex* z, int n_expc, dcomplex* z_expc) {
	printf("[INFO] %s applied to '%s'\n", solvername, funcname);
	int maxit = 30;
	double toler = TOLER;
	int verb = 0; //Verbosity level (0=quiet, 1=show iterations, 2=show also substeps).
	int i, conv[n];
	for (i = 0; i < n; i++) {//Loop on the roots to be tested.
		conv[i] = find_root(f, dom, &z[i], maxit, toler, verb);
		if (conv[i]) {
			printf("[INFO] %s converged to z[%d]=%g%+gi.\n", solvername, i, creal(z[i]), cimag(z[i]));
		}
		else {
			printf(STR_WARN "%s did not converged (initial z[%d]=%g%+gi).\n", solvername, i, creal(z[i]), cimag(z[i]));
		}
	}
	filter_cvector(&n, z, conv); //Removes the non-converged roots which are necessarily ignored for tests since we already know they have failed.
	print_cvector(n, z, "Found roots");
	return assert_cvector_subset(n, z, n_expc, z_expc, toler); //In principle, all roots in "z" must be actual roots of f(z)=0.
}

/**
 * Make all the test for the given root-finder.
 */
int test_find_root_method(char* solvername, int (*find_root)(dcomplex (*f)(dcomplex), Domain*, dcomplex*, int, double, int)) {
	print_test(solvername);
	int nf = 0;
	Domain* dom = (Domain*)calloc(1, sizeof(Domain));
	char* args[] = {"xrange=-7:+7", "yrange=-7:+7"};
	parse_domain(dom, sizeof(args)/sizeof(args[0]), args);
	
	char funcname1[] = "Simple polynomial function";
	dcomplex f1(dcomplex z) {
		return z*(z - 2. + 3.*I)*(z + 2. - 3.*I)*(z + 5);
	};
	dcomplex z1[] = {1-I, -5-5*I, 1, 6+3*I, 4, -6};
	dcomplex z1_expc[] = {0, 2-3*I, -2+3*I, -5};
	
	char funcname2[] = "Function with simple poles surrounding zero";
	dcomplex f2(dcomplex z) {
		return (z - 2.)/(z*(z - 2 + 2*I)*(z - 2 - 2*I)*(z - 4));
	};
	dcomplex z2[] = {1-I, -5-4*I, 1.5, 6+4*I, 4, -6};
	dcomplex z2_expc[] = {2};
	
	char funcname3[] = "Nearly constant function";
	dcomplex f3(dcomplex z) {
		return ctan(z*M_PI);
	};
	dcomplex z3[] = {2.5-I, -2+2*I, 1.8-1.5*I, 3+4*I, -4+3*I, 5+0.1*I}; //NB: The tangent is nearly constant far from the real axis.
	dcomplex z3_expc[] = {-7, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 7};
	
	nf += test_find_root_method_from_input(solvername, find_root, funcname1, f1, dom, sizeof(z1)/sizeof(z1[0]), z1, sizeof(z1_expc)/sizeof(z1_expc[0]), z1_expc);
	nf += test_find_root_method_from_input(solvername, find_root, funcname2, f2, dom, sizeof(z2)/sizeof(z2[0]), z2, sizeof(z2_expc)/sizeof(z2_expc[0]), z2_expc);
	nf += test_find_root_method_from_input(solvername, find_root, funcname3, f3, dom, sizeof(z3)/sizeof(z3[0]), z3, sizeof(z3_expc)/sizeof(z3_expc[0]), z3_expc);
	
	if (nf == 0) {
		printf(STR_PASS "%s: All tests passed.\n", solvername);
	}
	free(dom);
	return nf;
}

/**
 * Test the given solver "solve" with the given function "f(z)" for the given initial values "z" (array of size "n"),
 * and checks that all the resulting roots belong to the "z_expc" array. Returns 1 on failure, otherwise zero.
 */
int test_solve_method_from_input(char* solvername, void (*solve)(dcomplex (*f)(dcomplex), Domain*, int*, dcomplex*, int, double, int), char* funcname, dcomplex (*f)(dcomplex), Domain* dom, int n, dcomplex* z, int n_expc, dcomplex* z_expc) {
	printf("[INFO] %s applied to '%s'\n", solvername, funcname);
	int maxit = 30;
	double toler = TOLER;
	int verb = 0; //Verbosity level (0=quiet, 1=verbose, 2=debug).
	solve(f, dom, &n, z, maxit, toler, verb);
	print_cvector(n, z, "Found roots");
	return assert_cvector_subset(n, z, n_expc, z_expc, toler); //In principle, all roots in "z" must be actual roots of f(z)=0.
}

/**
 * Make all the test for the given root solver.
 */
int test_solve_method(char* solvername, void (*solve)(dcomplex (*f)(dcomplex), Domain*, int*, dcomplex*, int, double, int)) {
	print_test(solvername);
	int nf = 0;
	
	char funcname1[] = "Simple polynomial function";
	dcomplex f1(dcomplex z) {
		return (z - 2. + 3.*I)*(z + 2. - 3.*I)*(z + 5);
	};
	Domain dom1 = { .xmin = -10., .xmax = +10., .ymin = -10., .ymax = +10. };  //Search domain.
	int n1 = 5;  //Number of sought roots.
	dcomplex z1[n1];
	set_halton_points(&dom1, n1, z1);
	dcomplex z1_expc[] = {2-3*I, -2+3*I, -5};
	
	char funcname2[] = "Function with simple poles surrounding zero";
	dcomplex f2(dcomplex z) {
		return (z - 2.)/(z*(z - 2 + 2*I)*(z - 2 - 2*I)*(z - 4));
	};
	Domain dom2 = { .xmin = -1., .xmax = +5., .ymin = -3., .ymax = +3. };  //Search domain.
	int n2 = 5;  //Number of sought roots.
	dcomplex z2[n2];
	set_halton_points(&dom2, n2, z2);
	dcomplex z2_expc[] = {2};
	
	char funcname3[] = "Sine-like function";
	dcomplex f3(dcomplex z) {
		return ccos(z*M_PI/2);
	};
	Domain dom3 = { .xmin = -5.5, .xmax = +5.5, .ymin = -5.5, .ymax = +5.5 };  //Search domain.
	int n3 = 8;  //Number of sought roots.
	dcomplex z3[n3];
	set_halton_points(&dom3, n3, z3);
	dcomplex z3_expc[] = {-5, -3, -1, 1, 3, 5};
	//NB: Other tests are welcome...
	
	nf += test_solve_method_from_input(solvername, solve, funcname1, f1, &dom1, n1, z1, sizeof(z1_expc)/sizeof(z1_expc[0]), z1_expc);
	nf += test_solve_method_from_input(solvername, solve, funcname2, f2, &dom2, n2, z2, sizeof(z2_expc)/sizeof(z2_expc[0]), z2_expc);
	nf += test_solve_method_from_input(solvername, solve, funcname3, f3, &dom3, n3, z3, sizeof(z3_expc)/sizeof(z3_expc[0]), z3_expc);
	
	if (nf == 0) {
		printf(STR_PASS "%s: All tests passed.\n", solvername);
	}
	return nf;
}

/****************************************
 * TEST THE NONLINEAR FITTING ALGORITHM
 ***************************************/
/**
 * Test the Gauss-Newton fitting algorithm with the given data and compares with the expected result.
 */
int test_find_fit_gna_from_input(char* modelname, double (*f)(int i, double* p), int ns, double* ys, int np, double* param, int np_expc, double* param_expc) {
	printf("[TEST] Gauss-Newton fitting of '%s'...\n", modelname);
	int maxit = 30; //Default maximum number of Gauss-Newton iterations.
	double err_expc = 1.e-7;  //Expected relative error in the final fitting parameters.
	int verb = 0;  //Default verbosity level for tests.
	int conv = find_fit_gna(f, ns, ys, np, param, maxit, verb);
	if (conv == 0) {//If the method did not converged.
		printf(STR_WARN "The Gauss-Newton algorithm did not converged.\n");
		return 0; //In this case, we cannot say whether the algorithm worked or not, so returns OK flag.
	}
	print_vector(np, param, "Params");
	return assert_vector_close(np, param, np_expc, param_expc, err_expc);
}
 
/**
 * Make all the test of the Gauss-Newton fitting algorithm.
 */
int test_find_fit_gna() {
	print_test("find_fit_gna()");
	int nf = 0;
	
	double x_regular(int i, int ns) {//Local function to return the abscissa of the i-th point in the interval [0,1].
		double xmin = 0., xmax = 1., h = (xmax - xmin)/(ns - 1); //Interval of abscissa for all the 1D functions f(x).
		return xmin + i*h;
	}
	
	//Test #1: Linear model function.
	char modelname1[] = "Straight line";
	double ys1[] = {0.19, 0.77, 2.17, 2.81, 3.84, 4.74, 6.03, 6.84, 7.94, 9.12};
	double param1[] = {0., 0.}; //Initial parameters (a0, a1), with f(x) = a0 + a1*x.
	double param1_expc[] = {-0.040545454545454545, 8.97109090909090909};
	int ns1 = sizeof(ys1)/sizeof(ys1[0]);
	int np1 = sizeof(param1)/sizeof(param1[0]);
	double f1(int i, double* p) {//Model function.
		double x = x_regular(i, ns1);
		return p[0] + p[1]*x;
	}
	
	//Test #2: Sine wave model.
	char modelname2[] = "Sine wave";
	double ys2[] = {0.19, -0.018, 0.594, 0.418, 0.607, 0.623, 0.989, 0.837, 0.924, 1.049, 0.655, 0.837, 0.469, 0.168, 0.394,
		0.087, -0.064, -0.469, -0.812, -1.003, -0.881, -1.21, -1.115, -0.836, -0.949, -0.57, -0.71, -0.364, -0.204, -0.199};
	double param2[] = {1., 1., 1.}; //Initial amplitude "A", wave number "k", phase "phi".
	double param2_expc[] = {0.9749225716280133, 1.0005303292361316, 0.05253172519157753};
	int ns2 = sizeof(ys2)/sizeof(ys2[0]);
	int np2 = sizeof(param2)/sizeof(param2[0]);
	double f2(int i, double* p) {//Model function.
		double x = x_regular(i, ns2);
		return p[0]*sin(TWOPI*p[1]*x - p[2]);
	}
	
	//Test #3: Rational Pade [1/1] model function obtained by (3 + 2*x)/(1 + 6*x) plus a random perturbation of amplitude 0.1.
	char modelname3[] = "Pade [1/1]";
	double ys3[] = {3.044, 2.566, 2.367, 2.165, 2.001, 1.799, 1.653, 1.638, 1.58, 1.506, 1.44, 1.229, 1.355, 1.321, 1.149, 1.049, 1.104, 1.098, 1.069, 0.991, 
		1.02, 1.032, 0.955, 0.901, 0.849, 0.918, 0.851, 0.868, 0.737, 0.91, 0.862, 0.725, 0.877, 0.852, 0.728, 0.692, 0.809, 0.682, 0.671, 0.738};
	double param3[] = {2., 1., 10.};  //Initial Pade coefficients [a0, a1; b1], with b0=1 by convention.
	double param3_expc[] = {2.989809735787603, 2.0139951587983456, 5.945764148461985};
	int ns3 = sizeof(ys3)/sizeof(ys3[0]);
	int np3 = sizeof(param3)/sizeof(param3[0]);
	double f3(int i, double* p) {//Model function.
		double x = x_regular(i, ns3);
		return (p[0] + p[1]*x)/(1. + p[2]*x);
	}
	
	//Test #4: Continuous piecewise function based on the absolute value obtained by (1-2*abs(x-0.5)) plus a random perturbation of amplitude 0.2.
	char modelname4[] = "Absolute value";
	double ys4[] = {-0.009, -0.128, 0.077, 0.063, 0.148, 0.356, 0.445, 0.587, 0.513, 0.782, 0.88, 0.807, 0.667, 1.089,
		1.066, 0.819, 0.836, 0.701, 0.818, 0.637, 0.687, 0.444, 0.506, 0.536, 0.157, 0.26, 0.22, 0.154, 0.001, 0.03};
	double param4[] = {-1., 1., 0.3}; //Initial parameters of a*abs(x-c) + b.
	double param4_expc[] = {-2.1557357148487055, 1.0291511906217918, 0.4916965702716167};
	int ns4 = sizeof(ys4)/sizeof(ys4[0]);
	int np4 = sizeof(param4)/sizeof(param4[0]);
	double f4(int i, double* p) {//Model function.
		double x = x_regular(i, ns4);
		return p[0]*fabs(x - p[2]) + p[1];
	}
	
	//Test #5: Multi-dimensional spherically symmetric function based on the norm obtained by (1-2*norm(x-0.5)) plus a random perturbation of amplitude 0.2.
	char modelname5[] = "3D Norm function";
	double xs5[][3] = {{-0.043, -0.983, -0.306}, {-0.721, -0.639, 0.057}, {0.157, 0.521, -0.191}, {0.807, 0.952, 0.241}, {-0.804, 0.964, 0.5},
		{-0.733, -0.302, -0.635}, {0.296, -0.261, 0.33}, {-0.538, 0.117, 0.609}, {-0.938, -0.082, 0.064}, {0.081, -0.338, 0.152}, {0.187, -0.256, 0.083},
		{0.201, -0.75, -0.229}, {-0.609, -0.14, 0.829}, {-0.774, 0.995, 0.793}, {0.118, 0.664, -0.954}, {0.291, -0.802, -0.854}, {-0.789, 0.724, -0.74},
		{-0.003, -0.898, 0.403}, {-0.039, 0.427, -0.647}, {0.649, 0.639, -0.679}, {0.078, 0.951, 0.593}, {-0.62, 0.572, 0.475}, {0.21, 0.856, 0.571},
		{0.935, 0.549, -0.146}, {-0.538, -0.716, -0.871}, {0.012, -0.436, 0.538}, {-0.055, 0.47, 0.768}, {0.286, 0.159, 0.132}, {0.176, -0.625, -0.621},
		{-0.193, 0.998, 0.244}, {0.447, 0.215, 0.726}, {0.275, -0.592, 0.627}, {0.921, 0.488, 0.243}, {-0.893, 0.633, -0.358}, {0.277, 0.628, 0.819},
		{-0.263, 0.872, 0.699}, {-0.562, -0.143, 0.825}, {0.978, -0.567, -0.547}, {0.755, 0.994, -0.727}, {-0.46, -0.498, 0.}, {0.88, 0.046, 0.824},
		{-0.49, 0.168, 0.674}, {0.2, 0.616, -0.318}, {0.57, -0.003, 0.793}, {0.029, 0.749, -0.664}, {-0.646, 0.702, 0.849}, {0.816, -0.076, 0.452},
		{0.581, 0.692, 0.151}, {0.353, -0.927, -0.443}, {-0.853, 0.908, -0.663}};
	double ys5[] = {-1.156, -1.067, -0.173, -1.651, -1.661, -0.96, -0.018, -0.79, -0.785, 0.079, 0.157, -0.441, -1.2, -2.034, -1.39, -1.251, -1.44, 
		-0.877, -0.669, -1.426, -1.384, -0.922, -1.169, -1.277, -1.566, -0.303, -0.692, 0.292, -0.726, -0.902, -0.833, -0.655, -1.008, -1.123,
		-1.067, -1.216, -1.016, -1.532, -1.729, -0.348, -1.454, -0.582, -0.608, -0.968, -1.196, -1.471, -0.758, -0.788, -1.064, -1.641};
	double param5[] = {1., 1., 0.5, 0.5, 0.5};
	double param5_expc[] = {-1.8635625686127102, 0.8767736393808657, -0.007898780140093277, -0.014086372604185602, 0.0072762862506229814};
	int ns5 = sizeof(ys5)/sizeof(ys5[0]);
	int d = sizeof(xs5[0])/sizeof(xs5[0][0]); //Number of spatial dimensions, here 3.
	int np5 = sizeof(param5)/sizeof(param5[0]);
	if (ns5 != sizeof(xs5)/sizeof(xs5[0])) {
		fprintf(stderr, "[ERROR] Inconsistent number of samples, aborting...\n");
		return nf;
	}
	double f5(int i, double* p) {//Model function.
		return p[0]*distance(d, xs5[i], p+2) + p[1];
	}
	//NB: Other tests are welcome...
	
	nf += test_find_fit_gna_from_input(modelname1, f1, ns1, ys1, np1, param1, sizeof(param1_expc)/sizeof(param1_expc[0]), param1_expc);
	nf += test_find_fit_gna_from_input(modelname2, f2, ns2, ys2, np2, param2, sizeof(param2_expc)/sizeof(param2_expc[0]), param2_expc);
	nf += test_find_fit_gna_from_input(modelname3, f3, ns3, ys3, np3, param3, sizeof(param3_expc)/sizeof(param3_expc[0]), param3_expc);
	nf += test_find_fit_gna_from_input(modelname4, f4, ns4, ys4, np4, param4, sizeof(param4_expc)/sizeof(param4_expc[0]), param4_expc);
	nf += test_find_fit_gna_from_input(modelname5, f5, ns5, ys5, np5, param5, sizeof(param5_expc)/sizeof(param5_expc[0]), param5_expc);
	
	if (nf == 0) {
		printf(STR_PASS "find_fit_gna(): All tests passed.\n");
	}
	return nf;
}

/**
 * Main function
 */
int main(int argc, char** argv) {
	
	int nf = 0; //Current number of failed tests.
	
	nf += test_parametric_function();
	nf += test_eval_cpoly();
	nf += test_fit_poly();
	nf += test_fit_cpade();
	nf += test_find_root_poly();
	nf += test_find_root_method("find_root_secant()", find_root_secant);
	nf += test_find_root_method("find_root_muller()", find_root_muller);
	nf += test_solve_method("solve_from_guess()", solve_from_guess);
	nf += test_solve_method("solve_aberth()", solve_aberth_num);
	nf += test_find_fit_gna();
	
	printf("====== TEST RESULTS ======\n");
	if (nf == 0) {
		printf(STR_PASS "All tests are successful.\n");
	}
	else {
		printf(STR_FAIL "There are %d failed tests.\n", nf);
	}
	return 0;
}
