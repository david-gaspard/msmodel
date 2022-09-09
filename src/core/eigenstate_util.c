/****
 * @date Created on 2021-04-05 at 18:26:56 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to manipulate eigenstates.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "eigenstate_type.h"
#include "common_util.h"
#include "real_vector_util.h"
#include "complex_vector_util.h"
#include "medium_util.h"
#include "root_finder.h"

/**
 * Initialize an eigenstate with the given data.
 */
void init_eigenstate(Eigenstate* eigst, Medium* med) {
	eigst->eigv = (dcomplex*)calloc(med->n, sizeof(dcomplex)); //Allocate space for the components of the eigenvector.
	eigst->center = (double*)calloc(med->d, sizeof(double)); //Allocate space for the center.
}

/**
 * Free the content of the given eigenstate, not the pointer itself.
 */
void del_eigenstate(Eigenstate* eigst) {
	free(eigst->eigv);
	free(eigst->center);
}

/**
 * Computes the mass center of the given eigenstate weighted by the square modulus of the scattering amplitudes off the atoms.
 * The eigenstate center expected to be relatively close to the actual medium center.
 * The medium "med" is assumed to be initialized with the correct seed (to be checked by caller function).
 * The resulting center is written to the array "center", which is assumed to have "med->d" components.
 */
void center_eigenstate(Eigenstate* eigst, Medium* med, double* center) {//Private function.
	int c, d = med->d, n = med->n;
	for (c = 0; c < d; c++) {//First sets the center to zero.
		center[c] = 0.;
	}
	double re, im, a2, denom = 0.;
	int i;
	for (i = 0; i < n; i++) {//Loop on the atoms.
		re = creal(eigst->eigv[i]);
		im = cimag(eigst->eigv[i]);
		a2 = re*re + im*im;
		denom += a2;
		for (c = 0; c < d; c++) {//First sets the center to zero.
			center[c] += a2*med->pos[i*d + c];
		}
	}
	for (c = 0; c < d; c++) {//First sets the center to zero.
		center[c] /= denom;
	}
}

/**
 * Determines the list of the indices of the components of the given eigenstates "eigst" which are above the given tolerance "tol" in absolute value.
 * @param "eigst" An eigenstate, which is assumed to contains an eigenvector normlized to 1.
 * @param "tol" The tolerance on the components of the eigenstate, typically 1e-13.
 * @param "ns" The initial size of the array "idx". On exit, the final size of "idx".
 * @param "idx" Array of indices of valid eigenstate amplitudes above the tolerance "tol". Initially allocated with size "ns".
 * @returns "nexc" The number of excluded points which are below the tolerance "tol".
 */
int components_above_tol(Eigenstate* eigst, double tol, int* ns, int* idx) {//Private function.
	int i, size = 0;  //Number of components above the tolerance.
	for (i = 0; i < *ns; i++) {//Loop on the components of the eigenvector.
		if (cabs(eigst->eigv[i]) > tol) {
			idx[size] = i;
			size++;
		}
	}
	int nexc = *ns - size; //Number of excluded points.
	*ns = size;
	return nexc;
}

/**
 * Fits the given eigenstate "eigst" with the given medium "med" assuming the eigenvector resembles an exponential packet |a(x)| = A*exp(-locfac*|x-c|).
 * where "c" is the center of the eigenstate. The eigenvector is assumed to be normalized to one.
 */
void fit_eigenstate_exp(Eigenstate* eigst, Medium* med, int maxit, int verb) {
	//#1: Get all the valid components of the eigenvectors which are above the machine accuracy.
	int ns = med->n;
	int* idx = (int*)calloc(ns, sizeof(int));
	double tol = 1.e-13;  //Absolute lower tolerance on the modulus of the components of the eigenvector.
	int nexc = components_above_tol(eigst, tol, &ns, idx);
	if (verb >= 1 && nexc > 0)
		printf("[WARN] %d points are excluded from the fitting.\n", nexc);
	//#2: Change the representation of the data to logarithmic scale for more accurate results.
	double* ys = (double*)calloc(ns, sizeof(double));
	int i;
	for (i = 0; i < ns; i++) {//Fill the array of ordinates using only amplitudes above the machine accuracy.
		ys[i] = log(cabs(eigst->eigv[idx[i]]));
	}
	//#3: Defines the amplitude logarithm function which is fitted.
	double logamp(int i, double* p) {
		return p[0] - p[1]*distance(med->d, med->pos + idx[i]*med->d, &p[2]);
	};
	//#4: Initializes and roughly guess the parameters.
	int np = med->d+2; //Number of fitting parameters.
	double* param = (double*)calloc(np, sizeof(double));
	center_medium(med, eigst->center); //Guess the center using medium geometrical center.
	//center_eigenstate(eigst, med, eigst->center); //Old version using weighted average
	copy_vector(med->d, eigst->center, &param[2]);
	param[0] = 0.; param[1] = 0.1; //Nonzero values that does not really matter.
	//#5: Call the non-linear fitting function: find_fit();
	eigst->isfitted = find_fit_gna(logamp, ns, ys, np, param, maxit, verb);
	if (eigst->isfitted) {//If the fit is successful, then saves the results in the eigenstate.
		eigst->ampcst = param[0];
		eigst->locfac = param[1];
		copy_vector(med->d, &param[2], eigst->center);
	}
	free(idx);
	free(ys);
	free(param);
}

/**
 * Appends a "tikzpicture" environment to the given file "tikzfp".
 * The "tikzpicture" environment represents the absolute values of the scattering amplitudes of the given eigenstate "eigst"
 * on each atom as a function of the distance of these atoms from the mass center of the eigenstate. Gives to "eigst" the index "i".
 */
void plot_state_radial(Eigenstate* eigst, Medium* med, int i, FILE* tikzfp) {
	int n = med->n, d = med->d, i_hr = i+1; //One-based index for human readability.
	fprintf(tikzfp, "\\begin{tikzpicture}%%\n\\begin{axis}[%%\n\tname={st%d-rad},\n\tymode=log,\n\txmin=0, ymin=0,\n\ttitle={Amplitudes of the eigenstate $\\#%d$\\\\ at $k=(%g%+g\\,{\\rm i})\\,\\varsigma^{-1}$},\n\txlabel={Distance from center $r_i$~($\\varsigma$)},\n\tylabel={$|a_i|$~($\\varsigma^{%d}$)},\n]%%\n", i_hr, i_hr, creal(eigst->kroot), cimag(eigst->kroot), d-2);
	fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%st%d-rad\n", i_hr);
	double r, rmax = 0;
	int j;
	for (j = 0; j < n; j++) {//Loop on the components of the eigenvector, and the atoms of the medium.
		r = distance(d, eigst->center, med->pos + j*d); //Distance between the "center" and the "j"-th atom.
		fprintf(tikzfp, "\t(%.15g, %.15g)\n", r, cabs(eigst->eigv[j]));
		if (r > rmax) //Find the largest radius.
			rmax = r;
	}
	fprintf(tikzfp, "};\n");
	if (eigst->isfitted) {//If the fit has been successful, then show it graphically.
		fprintf(tikzfp, "\\addplot[blue, thick, mark=none, domain=0:%g, samples=10] {exp(%g - %g*x)};\n", rmax, eigst->ampcst, eigst->locfac);
	}
	else {
		printf("[WARN] The fitting of eigenstate #%d failed and could not be plotted.\n", i_hr);
	}
	fprintf(tikzfp, "\\end{axis}\n\\end{tikzpicture}%%\n");
}

/**
 * Appends a "tikzpicture" environment to the given file "tikzfp".
 * The "tikzpicture" environment represents the medium projected on the first two dimensions.
 * The radius of the disks depicts the absolute value of the scattering amplitudes of the given eigenstate "eigst".
 * Gives to "eigst" the zero-based index "i". This function assumes the eigenvectors are normalized to 1.
 */
void plot_state_proj(Eigenstate* eigst, Medium* med, int i, FILE* tikzfp) {
	int n = med->n, d = med->d, j, i_hr = i+1; //One-based index.
	fprintf(tikzfp, "\\begin{tikzpicture}%%\n\\begin{axis}[%%\n\tname={st%d-proj},\n", i_hr);
	if (d == 1) {//If d=1, then only shows the scattering amplitudes on each atom as vertical bars.
		fprintf(tikzfp, "\tymode=log,\n\tymin=0,\n\ttitle={Amplitudes of the eigenstate $\\#%d$\\\\ in a 1D medium},\n\txlabel={$x$~($\\varsigma$)},\n\tylabel={$|a_i|$~($\\varsigma^{%d}$)},\n]%%\n", i_hr, d-2);
		fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%st%d-proj\n\t", i_hr);
		for (j = 0; j < n; j++) {//Loop on the components of the eigenvector, and the atoms of the medium.
			fprintf(tikzfp, "(%.15g, %.15g) ", med->pos[j], cabs(eigst->eigv[j]));
		}
		fprintf(tikzfp, "\n};\n");
		//fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%center\n\t(%.15g, 0)\n};\n", eigst->center[0]);
	}
	else {//If d>=2, then use projected representation with disks.
		fprintf(tikzfp, "\ttitle={Amplitudes of the eigenstate $\\#%d$\\\\ in the %dD %s (projected view)},\n\txlabel={$x$~($\\varsigma$)},\n\tylabel={$y$~($\\varsigma$)},\n\taxis equal,\n\tdisabledatascaling,\n\tscatter src=explicit,\n\tscatter/@pre marker code/.code={%%\n\t\t\\pgfplotstransformcoordinatex{\\pgfplotspointmeta}%%\n\t\t\\scope[mark size=\\pgfplotsunitxlength*\\pgfmathresult]%%\n\t},\n]%%\n", i_hr, d, shape_to_string(med->shape));
		fprintf(tikzfp, "\\addplot[scatter, mark=*, only marks, red, fill opacity=0.5, draw opacity=0] coordinates {%%st%d-proj\n", i_hr);
		double r, ctr = 6.; //Contrast parameter, for instance, ctr=10.
		double rmin = 0.03, rmax = 0.40; //Minimum and maximum radius of disks assuming the maximum amplitude of an eigenvector component is 1, i.e., eigenvectors are normalized.
		for (j = 0; j < n; j++) {//Loop on the components of the eigenvector, and the atoms of the medium.
			r = rmin + (rmax - rmin)*asinh(sinh(ctr)*cabs(eigst->eigv[j]))/ctr; //Radius of the circles using asinh() function a safe logarithm against zero amplitude.
			fprintf(tikzfp, "\t(%.15g, %.15g) [%.15g]\n", med->pos[j*d], med->pos[j*d+1], r);
		}
		fprintf(tikzfp, "};\n");
		fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%center\n\t(%.15g, %.15g)\n};\n", eigst->center[0], eigst->center[1]);
	}
	fprintf(tikzfp, "\\end{axis}\n\\end{tikzpicture}%%\n");
}

/**
 * Appends a TikZ picture environment corresponding to the given eigenstate data contained
 * in "eigst" to the TikZ file "tikzfp". The plots are drawn on the same line. Gives to "eigst" the index "i".
 */
void export_tikz_eigenstate(Eigenstate* eigst, Medium* med, int i, FILE* tikzfp) {
	fprintf(tikzfp, "\\pgfplotsset{%%Plot eigenstate #%d\n\tevery node/.append style={font=\\footnotesize},\n\twidth=200pt,\n\ttitle style={align=center},\n\tevery node near coord/.style={font=\\tiny, inner sep=1.7pt},\n\tyticklabel style={rotate=90}\n}%%\n\\par%%\n", i+1);
	plot_state_proj(eigst, med, i, tikzfp); //Draw the eigenstate projected onto the two main dimensions.
	plot_state_radial(eigst, med, i, tikzfp); //Draw the eigenstate as a function of the radius from its center.
}

/**
 * Saves the given eigenstate stored in "eigst" to the file "fp". 
 * ses the index "i" to refer to "eigst" while saving.
 */
void save_eigenstate(Eigenstate* eigst, Medium* med, int i, FILE* fp) {
	fprintf(fp, "st%d_kroot=%.16lg%+.16lgi\n", i, creal(eigst->kroot), cimag(eigst->kroot));
	fprintf(fp, "st%d_mu=%.16lg%+.16lgi\n", i, creal(eigst->mu), cimag(eigst->mu));
	fprintf(fp, "st%d_eigv=", i);
	dcomplex z;
	int j;
	for (j = 0; j < med->n; j++) {//Loop on components of the eigenvector.
		z = eigst->eigv[j];
		fprintf(fp, "%.16lg%+.16lgi ", creal(z), cimag(z));
	}
	fprintf(fp, "\n");
}

/**
 * Parses the given eigenstate "eigst" from the given medium "med" and the "narg" arguments "args".
 * Uses the state index "i" to refer to "eigst" while parsing.
 * Note that this function assumes that the eigenstate has already been allocated with the appropriate init_*() function.
 */
void parse_eigenstate(Eigenstate* eigst, Medium* med, int i, int narg, char** args) {
	char key[20], *value;
	double re, im;
	sprintf(key, "st%d_kroot", i); //Parse the position of the eigenstate in the complex "k" plane.
	value = get_value(narg, args, key);
	if (sscanf(value, "%lg%lgi", &re, &im) != 2) {
		printf("[ERROR] Invalid k root '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	eigst->kroot = re + im*I; //Saves the eigenstate position in the complex "k" plane.
	sprintf(key, "st%d_mu", i); //Parse the eigenvalue.
	value = get_value(narg, args, key);
	if (sscanf(value, "%lg%lgi", &re, &im) != 2) {
		printf("[ERROR] Invalid eigenvalue '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	eigst->mu = re + im*I; //Saves the eigenvalue.
	sprintf(key, "st%d_eigv", i); //Parse the eigenvector.
	value = get_value(narg, args, key);
	if (!value[0]) {//If unexisting key "st<i>_eigv".
		printf("[ERROR] Data corruption: unexisting field '%s', aborting...\n", key);
		exit(EXIT_FAILURE);
	}
	int n = count_complex_data(value); //Count the complex numbers.
	if (n == med->n) {//If the number of elements is correct.
		parse_complex_data(n, eigst->eigv, value);
	}
	else {
		printf("[ERROR] Data corruption: found %d values in field '%s' but expected n=%d, aborting...\n", n, key, med->n);
		exit(EXIT_FAILURE);
	}
}
