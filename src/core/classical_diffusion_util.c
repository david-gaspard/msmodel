/****
 * @date Created on 2021-02-24 at 13:09:17 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the functions to determine the position of the classical transport pole in a d-ball
 * under the diffusion and hyperbolic diffusion approximations. Uses the multi-dimensional Newton-Raphson method.
 * Hyperbolic diffusion equations for (kappa, beta): 
 * (1) kappa^2/d + beta*(nsigma+beta) = 0
 * (2) (1/2d)*D_r I(d, kappa, R) + (V(d-1)/S(d))*(beta+nsigma)*I(d, kappa, R) = 0
 * Standard diffusion equations for (kappa, beta): 
 * (1) kappa^2/d + beta*nsigma = 0
 * (2) (1/2d)*D_r I(d, kappa, R) + (V(d-1)/S(d))*nsigma*I(d, kappa, R) = 0
 * where I(d,k,r) is the regular solution of the spherical Helmholtz equation, (nabla^2 + k^2)I(d,k,r) = 0, in dimension d.
 * One also uses the properties that: D_r I(d, k, r) = -2pi*r*I(d+2, k, r), and D_k I(d, k, r) = (k/2pi)*I(d-2, k, r).
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "green_bessel.h"
#include "real_vector_util.h"

#define PI     3.1415926535897932385   //Famous mathematical constant.
#define MEPS   1.0e-16                 //Approximate machine epsilon of double precision.
#define IMAX   50                      //Maximum number of terms allowed in the series.

/**
 * LAPACK's routine to compute the solution to a real system of linear equations.
 */
void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);

static int nvar = 2; //Number of sought variables of the diffusion equation (beta, kappa).

typedef struct DiffusionEquation_s {
	//Parameters:
	int d;          //Dimension.
	double radius;  //Radius of the d-ball medium, in units of "sp".
	double nsigma;  //Cross section of a single atom times the density (with unit density).
	//Unknowns:
	double beta;    //Variable dual to time (units of sp^-1).
	double kappa;   //Variable dual to radial coordinate (units of sp^-1).
	//Equation to be solved using Newton-Raphson iteration (computes the functions and the jacobian matrix):
	void (*equation)(struct DiffusionEquation_s* self, double* fun, double* jac);
} DiffusionEquation;

/**
 * Computes the real regular radial solution of the spherical Helmholtz equation, (nabla^2 + k^2) P(d, kr) = 0, in dimension "d".
 * This function is defined as P(d, kr) = HyperGeom_0F1(d/2, -(kr/2)^2), or in other words, Gamma(d/2)*(2/kr)^((d-2)/2)*J((d-2)/2, kr),
 * where J(nu, x) denotes the std Bessel function. The dimension "d" must be strictly positive (not zero).
 * The relative error less than 1e-13 up to k*r = 7, which is enough for the practical purposes, i.e., search for the fundamental mode.
 */
double helmholtz_series(int d, double k, double r) {
	double kr = k*r, x = -kr*kr/4;
	double t = 1., res = 1.;
	int i;
	for (i = 1; i < IMAX; i++) {
		t *= 2*x/(i*d);
		res += t;
		d += 2;
		if (fabs(t) < MEPS*fabs(res)) {
			//printf("[INFO] Break at step %d, t=%g, result=%g.\n", i, t, res);
			break;
		}
	}
	return res;
}

/**
 * Computes the "r"-derivative of the regular solution of the spherical Helmholtz equation in dimension "d".
 */
double dr_helmholtz(int d, double k, double r) {
	return -k*(k*r/d)*helmholtz_series(d+2, k, r);
}

/**
 * Computes the "k"-derivative of the regular solution of the spherical Helmholtz equation in dimension "d".
 */
double dk_helmholtz(int d, double k, double r) {
	return -r*(k*r/d)*helmholtz_series(d+2, k, r);
}

/**
 * Computes the derivative with respect to "k" of the "r"-derivative of the regular solution of the spherical Helmholtz equation in dimension "d".
 */
double dkdr_helmholtz(int d, double k, double r) {
	double kr = k*r;
	double P2 = helmholtz_series(d+2, k, r);
	double P4 = helmholtz_series(d+4, k, r);
	return (kr/d)*(kr*kr*P4/(d+2) - 2*P2);
}

/**
 * Returns the value of the functions to root-find and the Jacobian matrix for the standard diffusion equation.
 * The function values "f" are stored as {f(x,y), g(x,y)}, and the Jacobian matrix "j" in column-major format, for instance {D_x(f), D_x(g), D_y(f), D_y(g)}.
 */
void diffusion_equation_std(DiffusionEquation* deq, double* f, double* j) {
	double vds = uball_volume(deq->d-1)/uball_surface(deq->d);
	//Dispersion relation:
	f[0] = deq->kappa*deq->kappa/deq->d + deq->nsigma*deq->beta;
	//Boundary condition, (1/2d)*D_r P(d, kappa, R) + (V(d-1)/S(d))*nsigma*P(d, kappa, R) = 0 :
	f[1] = dr_helmholtz(deq->d, deq->kappa, deq->radius)/(2*deq->d) + vds*deq->nsigma*helmholtz_series(deq->d, deq->kappa, deq->radius);
	//Derivative of dipersion relation, f[0], with respect to "beta":
	j[0] = deq->nsigma;
	//Derivative of boundary condition, f[1], with respect to "beta":
	j[1] = 0;
	//Derivative of dispersion relation, f[0], with respect to "kappa":
	j[2] = 2*deq->kappa/deq->d;
	//Derivative of boundary condition, f[1], with respect to "kappa":
	j[3] = dkdr_helmholtz(deq->d, deq->kappa, deq->radius)/(2*deq->d) + vds*deq->nsigma*dk_helmholtz(deq->d, deq->kappa, deq->radius);
}

/**
 * Returns the value of the functions to root-find and the Jacobian matrix for the hyperbolic diffusion equation.
 * The function values "f" are stored as {f(x,y), g(x,y)}, and the Jacobian matrix "j" in column-major format, for instance {D_x(f), D_x(g), D_y(f), D_y(g)}.
 */
void diffusion_equation_hyp(DiffusionEquation* deq, double* f, double* j) {
	double vds = uball_volume(deq->d-1)/uball_surface(deq->d);
	double P = helmholtz_series(deq->d, deq->kappa, deq->radius);
	//Dispersion relation:
	f[0] = deq->kappa*deq->kappa/deq->d + (deq->beta + deq->nsigma)*deq->beta;
	//Boundary condition, (1/2d)*D_r P(d, kappa, R) + (V(d-1)/S(d))*nsigma*P(d, kappa, R) = 0 :
	f[1] = dr_helmholtz(deq->d, deq->kappa, deq->radius)/(2*deq->d) + vds*(deq->beta + deq->nsigma)*P;
	//Derivative of dipersion relation, f[0], with respect to "beta":
	j[0] = 2*deq->beta + deq->nsigma;
	//Derivative of boundary condition, f[1], with respect to "beta":
	j[1] = vds*P;
	//Derivative of dispersion relation, f[0], with respect to "kappa":
	j[2] = 2*deq->kappa/deq->d;
	//Derivative of boundary condition, f[1], with respect to "kappa":
	j[3] = dkdr_helmholtz(deq->d, deq->kappa, deq->radius)/(2*deq->d) + vds*(deq->beta + deq->nsigma)*dk_helmholtz(deq->d, deq->kappa, deq->radius);
}

/**
 * Returns the value of the functions to root-find and the Jacobian matrix for the diffusion equation (version system 1, "S1").
 * It is essentially the same as the standard diffusion approximation but with a difference in the boundary condition.
 * The function values "f" are stored as {f(x,y), g(x,y)}, and the Jacobian matrix "j" in column-major format, for instance {D_x(f), D_x(g), D_y(f), D_y(g)}.
 */
void diffusion_equation_sys1(DiffusionEquation* deq, double* f, double* j) {
	double sd2v = uball_surface(deq->d)/(2*uball_volume(deq->d-1));
	//Dispersion relation:
	f[0] = deq->kappa*deq->kappa/deq->d + deq->nsigma*deq->beta;
	//Boundary condition, D_r P(d, kappa, R) + (S(d)/2V(d-1))*nsigma*P(d, kappa, R) = 0 :
	f[1] = dr_helmholtz(deq->d, deq->kappa, deq->radius) + sd2v*deq->nsigma*helmholtz_series(deq->d, deq->kappa, deq->radius);
	//Derivative of dipersion relation, f[0], with respect to "beta":
	j[0] = deq->nsigma;
	//Derivative of boundary condition, f[1], with respect to "beta":
	j[1] = 0;
	//Derivative of dispersion relation, f[0], with respect to "kappa":
	j[2] = 2*deq->kappa/deq->d;
	//Derivative of boundary condition, f[1], with respect to "kappa":
	j[3] = dkdr_helmholtz(deq->d, deq->kappa, deq->radius) + sd2v*deq->nsigma*dk_helmholtz(deq->d, deq->kappa, deq->radius);
}

/**
 * Returns the value of the functions to root-find and the Jacobian matrix for the diffusion equation (version system 2, "S2").
 * It is essentially the same as the hyperbolic diffusion approximation but with some differences.
 * The function values "f" are stored as {f(x,y), g(x,y)}, and the Jacobian matrix "j" in column-major format, for instance {D_x(f), D_x(g), D_y(f), D_y(g)}.
 */
void diffusion_equation_sys2(DiffusionEquation* deq, double* f, double* j) {
	double sd2v = uball_surface(deq->d)/(2*uball_volume(deq->d-1));
	double P = helmholtz_series(deq->d, deq->kappa, deq->radius);
	//Dispersion relation:
	f[0] = deq->nsigma*deq->kappa*deq->kappa/deq->d + deq->beta*(deq->beta + deq->nsigma)*(deq->beta + deq->nsigma);
	//Boundary condition, D_r P(d, kappa, R) + (S(d)/2V(d-1))*nsigma*(1 + 2*beta/nsigma)*(1 + beta/nsigma)*P(d, kappa, R) = 0 :
	f[1] = deq->nsigma*dr_helmholtz(deq->d, deq->kappa, deq->radius) + sd2v*(2*deq->beta + deq->nsigma)*(deq->beta + deq->nsigma)*P;
	//Derivative of dipersion relation, f[0], with respect to "beta":
	j[0] = (3*deq->beta + deq->nsigma)*(deq->beta + deq->nsigma);
	//Derivative of boundary condition, f[1], with respect to "beta":
	j[1] = sd2v*(4*deq->beta + 3*deq->nsigma)*P;
	//Derivative of dispersion relation, f[0], with respect to "kappa":
	j[2] = 2*deq->nsigma*deq->kappa/deq->d;
	//Derivative of boundary condition, f[1], with respect to "kappa":
	j[3] = deq->nsigma*dkdr_helmholtz(deq->d, deq->kappa, deq->radius) + sd2v*(2*deq->beta + deq->nsigma)*(deq->beta + deq->nsigma)*dk_helmholtz(deq->d, deq->kappa, deq->radius);
}

/**
 * Solve the given system of diffusion equations for (beta, kappa) using the Newton-Raphson iteration in two dimensions.
 * The unknowns (beta, kappa) must be initialized with a reasonable guess.
 */
int newton_solve(DiffusionEquation* deq, int maxit, double toler, int verbose, const char* msg) {
	if (verbose >= 2) {
		printf("[INFO] %s | Initial values: beta=%.16g,\tkappa=%.16g \n", msg, deq->beta, deq->kappa);
	}
	int nrhs = 1, ipiv[nvar], info, s, conv = 0;  //By default, no convergence.
	double x[nvar], f[nvar], j[nvar*nvar];
	x[0] = deq->beta;
	x[1] = deq->kappa;
	for (s = 1; s <= maxit; s++) {//Newton-Raphson iteration for "kappa".
		deq->equation(deq, f, j);
		if (verbose >= 2) {
			printf("[INFO] %s | Step #%02d | f=%g,\t", msg, s, norm(nvar, f));
		}
		dgesv_(&nvar, &nrhs, j, &nvar, ipiv, f, &nvar, &info);
		add_vector(nvar, x, f, -1);
		deq->beta  = x[0]; //Saves the result.
		deq->kappa = x[1];
		if (verbose >= 2) {
			printf("beta=%.16g,\tkappa=%.16g,\tres=%g \n", deq->beta, deq->kappa, norm(nvar, f));
		}
		if (norm(nvar, f) <= toler*norm(nvar, x)) {
			conv = 1;
			break;
		}
	}
	if (verbose >= 1) {
		if (conv) {
			printf("[INFO] %s | Found root beta=%g sp^-1, kappa=%g sp^-1, after %d steps (final res=%g).\n", msg, deq->beta, deq->kappa, s, norm(nvar, f));
		}
		else {
			printf("[WARN] %s | No convergence after %d steps (final res=%g). The result is not meaningful.\n", msg, s, norm(nvar, f));
		}
	}
	return conv;
}

/**
 * Solve the classical transport equations under the diffusion approximation and the hyperbolic diffusion approximation.
 * Saves the values of "beta" and "kappa" associated with the diffusion equation, and "betahyp" and "kappahyp" associated with the hyperbolic diffusion equation.
 * The flags "conv" and "convhyp" are equal to 1 only if the Newton-Raphson iteration converged for the standard and hyperbolic equations, respectively.
 */
void diffusion_solve_v1(int d, double radius, double nsigma, double* beta, double* kappa, int* conv, double* betahyp, double* kappahyp, int* convhyp) {
	//STEP 0: Initialize the diffusion equations.
	DiffusionEquation* deq = (DiffusionEquation*)calloc(1, sizeof(DiffusionEquation));
	deq->d = d;
	deq->radius = radius;
	deq->nsigma = nsigma;
	//STEP 1: Initialize the unknowns (beta, kappa) with appropriate guess.
	deq->beta  = 0;
	deq->kappa = sqrt(2*deq->d)/deq->radius;
	//STEP 2: Solve the standard diffusion equation.
	int maxit = 20;
	double toler = 1e-13;
	int verbose = 1;
	deq->equation = diffusion_equation_std;
	*conv  = newton_solve(deq, maxit, toler, verbose, "Std");
	*beta  = deq->beta;   //Saves the solution of the standard diffusion.
	*kappa = deq->kappa;
	//STEP 3: Solve the hyperbolic diffusion equation.
	deq->equation = diffusion_equation_hyp;  //Change of equation, so that start from the previous guess (beta, kappa).
	*convhyp  = newton_solve(deq, maxit, toler, verbose, "Hyp");
	*betahyp  = deq->beta;   //Saves the solution of the hyperbolic diffusion.
	*kappahyp = deq->kappa;
	free(deq); //Frees the workspace.
}

/**
 * Writes the TikZ option to plot the solution of the diffusion equation in the string "str".
 * The output string "str" must be allocated with at least 600 characters.
 * Returns the final length of the string "str".
 * This old version uses the convention of the original draft notes of the thesis.
 */
int diffusion_result_plot_v1(int d, double radius, double nsigma, char* str) {
	double beta, kappa, betahyp, kappahyp, pos;
	int conv, convhyp;
	diffusion_solve_v1(d, radius, nsigma, &beta, &kappa, &conv, &betahyp, &kappahyp, &convhyp);
	char cmd1[300] = {0}, cmd2[300] = {0};
	if (conv) {
		printf("[INFO] Plotting the standard diffusion pole at beta=%g sp^-1...\n", beta);
		pos = -beta/2;  //The factor of 2 comes from the classical vs quantum correspondence.
		sprintf(cmd1, "\\draw[black, densely dotted, thick] (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[black!70,pos=0,above right,rotate=90]{DA};\n", pos, pos);
	}
	if (convhyp) {
		printf("[INFO] Plotting the hyperbolic diffusion pole at beta=%g sp^-1...\n", betahyp);
		pos = -betahyp/2;
		sprintf(cmd2, "\\draw[black, densely dashed, thick] (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[black!70,pos=0,below right,rotate=90]{HDA};\n", pos, pos);
	}
	return sprintf(str, "%s%s", cmd1, cmd2);
}

/**
 * Writes the TikZ option to plot the solution of the diffusion equation in the string "str".
 * The output string "str" must be allocated with at least 600 characters.
 * Returns the final length of the string "str".
 * This version uses the conventions of the final report of the thesis.
 */
int diffusion_result_plot_v2(int d, double radius, double nsigma, char* str) {
	double beta1, kappa1, beta2, kappa2, pos;
	int conv1, conv2;
	//STEP 1: Initialize the diffusion equations.
	DiffusionEquation* deq = (DiffusionEquation*)calloc(1, sizeof(DiffusionEquation));
	deq->d = d;
	deq->radius = radius;
	deq->nsigma = nsigma;
	deq->kappa = sqrt(2*deq->d)/deq->radius;
	deq->beta  = -deq->kappa*deq->kappa/(deq->d*deq->nsigma);
	//STEP 2: Solve the diffusion equation version "S1".
	int maxit = 20;
	double toler = 1e-13;
	int verbose = 1; //Enable verbosity.
	deq->equation = diffusion_equation_sys1;
	conv1  = newton_solve(deq, maxit, toler, verbose, "S1");
	beta1  = deq->beta;   //Saves the solution of the standard diffusion.
	kappa1 = deq->kappa;
	//STEP 3: Solve the diffusion equation version "S2".
	deq->equation = diffusion_equation_sys2;  //Change of equation, so that start from the previous guess (beta, kappa).
	conv2  = newton_solve(deq, maxit, toler, verbose, "S2");
	beta2  = deq->beta;   //Saves the solution of the hyperbolic diffusion.
	kappa2 = deq->kappa;
	free(deq); //Frees the workspace.
	//STEP 4: Plot the results in the string "str".
	char cmd1[300] = {0}, cmd2[300] = {0};
	if (conv1) {
		printf("[INFO] Plotting the S1 diffusion pole at beta=%.12lg sp^-1, kappa=%.12lg sp^-1...\n", beta1, kappa1);
		pos = -beta1/2;  //The factor of 2 comes from the classical vs quantum correspondence.
		sprintf(cmd1, "\\draw[black, densely dotted, thick] (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[black!60,pos=0,above right,rotate=90]{S1};\n", pos, pos);
	}
	if (conv2) {
		printf("[INFO] Plotting the S2 diffusion pole at beta=%.12lg sp^-1, kappa=%.12lg sp^-1...\n", beta2, kappa2);
		pos = -beta2/2;
		sprintf(cmd2, "\\draw[black, densely dashed, thick] (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[black!60,pos=0,below right,rotate=90]{S2};\n", pos, pos);
	}
	return sprintf(str, "%s%s", cmd1, cmd2);
}
