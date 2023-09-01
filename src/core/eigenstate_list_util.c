/****
 * @date Created on 2021-03-29 at 19:42:54 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the eigenstate utilities.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "common_util.h"
#include "real_vector_util.h"
#include "complex_vector_util.h"
#include "domain_util.h"
#include "root_finder.h"
#include "medium_util.h"
#include "eigenstate_util.h"
#include "eigenstate_list_type.h"

/**
 * Initializes the eigenstate object with least number of parameters.
 * Note that the eigenvectors are allocated later while parsing or just before computations.
 */
void init_eigenstate_list(EigenstateList* eigls, Medium* med, uint64_t seed, EigenstateMethod method, int ntarget, Domain* dom, int maxit, double toler, int verb) {
	if (ntarget <= 0 || ntarget > 1e4) {
		printf("[ERROR] Invalid number of root-finding targets, ntarget=%d, aborting...\n", ntarget);
		exit(EXIT_FAILURE);
	}
	//if (nexpsup < 0 || nexpsup > 1e2) {
	//	printf("[ERROR] Invalid number of exponential suppression parameters, nexpsup=%d, aborting...\n", nexpsup);
	//	exit(EXIT_FAILURE);
	//}
	if (maxit <= 2 || maxit > 1e4) {
		printf("[ERROR] Invalid number of root-finding iterations, maxit=%d, aborting...\n", maxit);
		exit(EXIT_FAILURE);
	}
	if (toler <= 0. || toler > 0.5) {
		printf("[ERROR] Invalid tolerance, toler=%g, aborting...\n", toler);
		exit(EXIT_FAILURE);
	}
	eigls->med = med;
	eigls->seed = seed;
	fill_medium(med, seed); //Fill medium here, assuming only one geometry is expected.
	eigls->method = method;
	eigls->ntarget = ntarget;
	eigls->dom = (Domain*)calloc(1, sizeof(Domain));
	copy_domain(dom, eigls->dom);
	eigls->maxit = maxit;
	eigls->toler = toler;
	eigls->verb = verb;
}

/**
 * Allocate space for the eigenstates.
 */
void alloc_eigenstate_list(EigenstateList* eigls, int nstate) {
	if (nstate < 0 || nstate > eigls->ntarget) {
		printf("[ERROR] Invalid number of states, nstate=%d, expected less than %d states, aborting...\n", nstate, eigls->ntarget);
		exit(EXIT_FAILURE);
	}
	eigls->nstate = nstate;
	eigls->eigtab = (Eigenstate**)calloc(nstate, sizeof(Eigenstate*));
	int i;
	for (i = 0; i < nstate; i++) {//Loop on all the states to initialize them individually.
		eigls->eigtab[i] = (Eigenstate*)calloc(1, sizeof(Eigenstate));
		init_eigenstate(eigls->eigtab[i], eigls->med);
	}
}

/**
 * Free the content of eigenstate, not the pointer itself.
 */
void del_eigenstate_list(EigenstateList* eigls) {
	int i;
	for (i = 0; i < eigls->nstate; i++) {//Loop on all the states to free each of them.
		del_eigenstate(eigls->eigtab[i]);
		free(eigls->eigtab[i]);
	}
	free(eigls->eigtab);
	free(eigls->dom);
}

/**
 * Converts the given string to the corresponding eigenstate method.
 * Returns 1 if the string is not recognized. 
 */
int string_to_eigmethod(const char* str, EigenstateMethod* method) {
	if (strcmp(str, "determinant") == 0) {
		*method = determinant;
	}
	else if (strcmp(str, "mineigval") == 0) {
		*method = mineigval;
	}
	else if (strcmp(str, "invtraceinv") == 0) {
		*method = invtraceinv;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Converts the given eigenstate method to the corresponding string.
 * Returns the string "none" if the method is not recognized. 
 */
char* eigmethod_to_string(EigenstateMethod method) {
	switch (method) {
		case determinant:
			return "determinant";
		case mineigval:
			return "mineigval";
		case invtraceinv:
			return "invtraceinv";
		default:
			return "none";
	}
}

/**
 * Check the eigenstates contained in the given eigenstates list "eigls",
 * in particular if there is no duplicate states and if their eigenvalues are small enough.
 * This function assumes that the eigenstates exist in memory.
 */
void check_eigenstate_list(EigenstateList* eigls) {
	int issuspect, nsuspect = 0;  //Current number of suspect states.
	double delta, toldk, tolmu = 1e-6; //Tolerance over the "mu" eigenvalue in absolute value.
	dcomplex mui, krooti, krootj;
	int i, j;
	for (i = 0; i < eigls->nstate; i++) {//Loop again on the states to inform the user about their convergence.
		issuspect = 0;
		mui = eigls->eigtab[i]->mu;
		if (cabs(mui) > tolmu) {//Warns the user in case of non-convergence, but expected to be exceptionally rare.
			issuspect = 1;
			printf("[WARN] Eigenstate %d/%d is probably unreliable: abs(mu)=%g (tolmu=%g).\n",
				(i+1), eigls->nstate, cabs(mui), tolmu);
		}
		krooti = eigls->eigtab[i]->kroot;
		toldk = tolmu*cabs(krooti);
		for (j = 0; j < i; j++) {//Loop on the previous state for any posible duplicate.
			krootj = eigls->eigtab[j]->kroot;
			delta = cabs(krooti - krootj);
			if (delta < toldk) {
				issuspect = 1;
				printf("[WARN] Eigenstates %d/%d and %d/%d are probably identical, k%d=%.15g%+.15gi, k%d=%.15g%+.15gi, delta=%g.\n",
					(i+1), eigls->nstate, (j+1), eigls->nstate, i, creal(krooti), cimag(krooti), j, creal(krootj), cimag(krootj), delta);
			}
		}
		nsuspect += issuspect;
	}
	if (nsuspect != 0) {
		printf("[WARN] There are %d probable spurious states.\n", nsuspect);
	}
	else {
		printf("[INFO] All states are reliable and distinct.\n");
	}
}

/**
 * Fit all the eigenstates in the given eigenstate list "eigls".
 * This function assumes that the medium "eigls->med" is properly filled with the correct seed "eigls->seed".
 */
void fit_eigenstate_list(EigenstateList* eigls, int maxit, int verb) {//Private function.
	int i, nfailed = 0;
	int* ifailed = (int*)calloc(eigls->nstate, sizeof(int));  //Indices of failed states.
	for (i = 0; i < eigls->nstate; i++) {//Loop on the states to plot them.
		if (verb >= 1) {
			printf("[INFO] Fitting model to the eigenstate %3d/%3d ...\n", (i+1), eigls->nstate); fflush(stdout);
		}
		fit_eigenstate_exp(eigls->eigtab[i], eigls->med, maxit, verb); //Fit the eigenstates assuming the medium is filled with the correct seed.
		eigls->nfitted += eigls->eigtab[i]->isfitted; //Saves the number of successfully fitted states.
		if (eigls->eigtab[i]->isfitted != 1) {//If the fitting has failed, then save the index of the concerned state.
			ifailed[nfailed] = i;
			nfailed++; //Counts the number of failed fitttings.
		}
	}
	if (nfailed != 0) {
		printf("[WARN] Failed to fit %d eigenstates: ", nfailed);
		for (i = 0; i < nfailed; i++) {//Display the indices of eigenstates that have been failed to fit.
			printf("#%d ", (ifailed[i]+1));
		}
		printf("\n");
	}
	else {
		printf("[INFO] All states have been successfully fitted.\n");
	}
	fflush(stdout);
	free(ifailed);
}

/**
 * Plot the found positions of the roots to det(M(k))=0 in the complex "k" plane in the given file "tikzfp".
 */
void plot_state_kroots(EigenstateList* eigls, FILE* tikzfp) {//Private function.
	fprintf(tikzfp, "\\begin{axis}[%%\n\tname={roots},\n\tymin=%g, ymax=%g,\n\ttitle={(a) Resonances},\n\txlabel={$-\\Im k~(\\varsigma^{-1})$},\n\tylabel={$\\Re k~(\\varsigma^{-1})$},\n\txticklabel pos=right,\n]%%\n", eigls->dom->xmin, eigls->dom->xmax);
	fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%Resonances\n"); //Add this to show indices: "nodes near coords, point meta=explicit symbolic".
	dcomplex kroot;
	int i;
	for (i = 0; i < eigls->nstate; i++) {//Loop on the eigenstates.
		kroot = eigls->eigtab[i]->kroot;
		fprintf(tikzfp, "\t(%.15g, %.15g) [%d]\n", -cimag(kroot), creal(kroot), (i+1));
	}
	fprintf(tikzfp, "};\n\\end{axis}%%\n");
}

/**
 * Plot the localization factors for each eigenstate in the given file "tikzfp".
 */
void plot_state_locfac(EigenstateList* eigls, FILE* tikzfp) {//Private function.
	int i, size = eigls->nstate; //The number of points is initially the number of eigenstates.
	int ns = 0; //Actual number of sample points to be fitted and plotted, i.e., logical size of the arrays "xs" and "ys".
	double xs[size], ys[size];
	for (i = 0; i < size; i++) {//Loop on the states to initialize the data points
		if (eigls->eigtab[i]->isfitted) {//If the state has been successfully fitted, then the localization factor is meaningful.
			xs[ns] = -cimag(eigls->eigtab[i]->kroot); //Abscissa: -Im(k).
			ys[ns] = eigls->eigtab[i]->locfac; //Ordinate: Lambda.
			ns++; //Increment the actual number of points.
		}
	}
	int d = 1; //Degree of the polynomial fit.
	double poly[d+1];
	fit_poly(ns, xs, ys, d, poly); //Simple linear fitting of the localization factor.
	fprintf(tikzfp, "\\begin{axis}[%%\n\tname={locfac},\n\tat={(roots.south west)},\n\ttitle={(b) Localization factor},\n\tylabel={$\\Lambda~(\\varsigma^{-1})$},\n\textra y ticks={0},\n\textra y tick style={grid=major},\n\tlegend pos=south east,\n\tlegend cell align=left,\n\txticklabels=none,\n]%%\n");
	fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%Localization factor\n");
	for (i = 0; i < ns; i++) {//Loop on the eigenstates to plot them.
		fprintf(tikzfp, "\t(%.15g, %.15g) [%d]\n", xs[i], ys[i], (i+1));
	}
	fprintf(tikzfp, "}; \\addlegendentry{data points $(-k_{\\rm i},\\Lambda)$}\n");
	fprintf(tikzfp, "\\addplot[blue, thick, domain=%g:%g, samples=5] {%.15g %+.15g*x}; \\addlegendentry{$%.4f(-k_{\\rm i})%+.4f$}\n\\end{axis}%%\n",
		-eigls->dom->ymax, -eigls->dom->ymin, poly[0], poly[1], poly[1], poly[0]);
}

/**
 * Plot the distances between the eigenstate center and the medium center for each eigenstate in the given file "tikzfp".
 */
void plot_state_center(EigenstateList* eigls, FILE* tikzfp) {//Private function.
	fprintf(tikzfp, "\\begin{axis}[%%\n\tname={center},\n\tat={(locfac.south west)},\n\tymin=0, ymax=2,\n\trestrict y to domain=0:2,\n\ttitle={(c) Centers of states},\n\txlabel={$-\\Im k~(\\varsigma^{-1})$},\n\tylabel={$\\|\\vec{c}\\|/R$},\n]%%\n");
	fprintf(tikzfp, "\\addplot[mark=*, only marks, mark options={black, scale=0.3}] coordinates {%%Centers of states\n");
	int i, d = eigls->med->d; //Number of spatial dimensions.
	dcomplex kroot;
	double center[d], dist, radius = radius_of_ball(eigls->med); //Compute the radius of the medium as if it is a d-ball.
	center_medium(eigls->med, center); //Compute the geometrical center of the medium.
	for (i = 0; i < eigls->nstate; i++) {//Loop on the eigenstates.
		kroot = eigls->eigtab[i]->kroot;
		dist = distance(d, center, eigls->eigtab[i]->center); //Compute the distance between the eigenstate center and the medium center.
		fprintf(tikzfp, "\t(%.15g, %.15g) [%d]\n", -cimag(kroot), dist/radius, (i+1));
	}
	fprintf(tikzfp, "};\n\\end{axis}%%\n");
}

/**
 * Draw to the file "tikzfp" several plots summarizing the properties of the found eigenstates such as: (a) the positions of the resonances, (b) the localization factor, and (c) the relative distance between the eigenstate center and the medium center.
 */
void plot_summary(EigenstateList* eigls, FILE* tikzfp) {//Private function.
	fprintf(tikzfp, "\\begin{center}%%\n%dD %s with $N=%d$ atoms and seed %lu.\\\\ Found %d/%d complex roots in %d iterations with method ``%s'',\\\\ of which %d successfully fitted states.\\\\ Total computation time is %g s (avg %g root/s).\n\\end{center}%%\n",
		eigls->med->d, shape_to_string(eigls->med->shape), eigls->med->n, eigls->seed, eigls->nstate, eigls->ntarget, eigls->maxit, eigmethod_to_string(eigls->method), eigls->nfitted, eigls->realtime, eigls->nstate/eigls->realtime);
	fprintf(tikzfp, "\\begin{tikzpicture}[%%\n\t/pgfplots/every axis/.append style={xmin=%g, xmax=%g, width=300pt, height=170pt, anchor=north west, yshift=-10pt},\n\t/pgfplots/every axis title/.style={below right, at={(0,1)}, fill=white, fill opacity=0.7, text opacity=1, inner sep=1.2pt},\n]%%\n", -eigls->dom->ymax, -eigls->dom->ymin);
	plot_state_kroots(eigls, tikzfp);
	plot_state_locfac(eigls, tikzfp);
	plot_state_center(eigls, tikzfp);
	fprintf(tikzfp, "\\end{tikzpicture}%%\n");
}

/**
 * Appends several TikZ picture environments corresponding to all tha data of the given eigenstate list object "eigls"
 * to the file of given name "fname". Note that the filename "fname" should not have any extension.
 */
void export_tikz_eigenstate_list(EigenstateList* eigls, const char* fname) {
	int maxit = 100; //Maximum number of Gauss-Newton fitting iterations. The convergence may be slow when the exponential decrease is not significant.
	int verb = 0;    //Verbosity level of the Gauss-Newton algorithm (0=quiet, 1=verbose). Can be used to debug.
	fit_eigenstate_list(eigls, maxit, verb); //First fit all the eigenstates.
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	fprintf(tikzfp, "\\pgfplotsset{%%\n\tevery axis/.append style={%%\n\t\tfont=\\footnotesize,\n\t\twidth=180pt,\n\t\ttitle style={align=center},\n\t\tevery node near coord/.style={font=\\tiny, inner sep=1.7pt},\n\t\tyticklabel style={rotate=90}\n\t}\n}%%\n");
	plot_summary(eigls, tikzfp); //Draw the summary plots (including localization factors).
	int i, iplot, nplot = 5; //Upper limit on the number of eigenstates to be plotted.
	if (eigls->nstate < nplot) nplot = eigls->nstate; //If the number of state is lower than "nplot", then plot all the eigenstates.
	for (iplot = 0; iplot < nplot; iplot++) {//Loop on a few states to plot them.
		i = (nplot != 1) ? ((eigls->nstate - 1)*iplot)/(nplot - 1) : 0; //Actual zero-based index of the eigenstate to plot.
		fprintf(tikzfp, "\\par%%\n");
		plot_state_proj(eigls->eigtab[i], eigls->med, i, tikzfp);
		plot_state_radial(eigls->eigtab[i], eigls->med, i, tikzfp);
	}
	fclose(tikzfp);
}

/**
 * Saves all the data contained in the eigenstate list "eigls" to the file of given data name "fname".
 */
void save_eigenstate_list(EigenstateList* eigls, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[eigenstate_list]\n");  //Seed used for the random medium.
	fprintf(fp, "seed=%lu\n", eigls->seed);
	fprintf(fp, "method=%s\n", eigmethod_to_string(eigls->method));
	fprintf(fp, "ntarget=%d\n", eigls->ntarget);
	fprintf(fp, "nstate=%d\n", eigls->nstate);
	save_domain(eigls->dom, fp);
	fprintf(fp, "maxit=%d\n", eigls->maxit);
	fprintf(fp, "toler=%g\n", eigls->toler);
	fprintf(fp, "verb=%d\n", eigls->verb);
	fprintf(fp, "realtime=%g\n", eigls->realtime); //Total physical duration of the computation (in seconds).
	int i;
	for (i = 0; i < eigls->nstate; i++) {//Loop on the unique states to save them.
		save_eigenstate(eigls->eigtab[i], eigls->med, i, fp);
	}
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Parses the eigenstate list "eigls" from the "narg" arguments contained in "args".
 * The "eigls" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_eigenstate_list(EigenstateList* eigls, Medium* med, int narg, char** args) {
	int nargmin = 4; //Minimum number of arguments (seed/maxit/toler/ntarget).
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	uint64_t seed = 0;  //Default values.
	EigenstateMethod method;
	int nstate, ntarget = 0, maxit = 0, verb = 1;
	double toler = 0.;
	char* value;
	value = get_value(narg, args, "seed");
	if (sscanf(value, "%lu", &seed) != 1) {
		printf("[ERROR] Invalid seed '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "method");
	if (string_to_eigmethod(value, &method)) {
		printf("[ERROR] Invalid eigenstate method '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "ntarget");
	if (sscanf(value, "%d", &ntarget) != 1) {
		printf("[ERROR] Invalid number of targets '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "maxit");
	if (sscanf(value, "%d", &maxit) != 1) {
		printf("[ERROR] Invalid number of iterations '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "toler");
	if (sscanf(value, "%lg", &toler) != 1) {
		printf("[ERROR] Invalid tolerance '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	sscanf(get_value(narg, args, "verb"), "%d", &verb);
	Domain dom; //Parse the complex domain.
	if (parse_domain(&dom, narg, args)) {//If the domain parsing failed
		printf("[ERROR] Invalid domain, aborting...\n");
		exit(EXIT_FAILURE);
	}
	init_eigenstate_list(eigls, med, seed, method, ntarget, &dom, maxit, toler, verb);
	value = get_value(narg, args, "nstate");
	if (sscanf(value, "%d", &nstate) == 1) {//If existing field "nstate", then also parse the eigenstates.
		alloc_eigenstate_list(eigls, nstate); //Allocate space for all the eigenstates.
		int i;
		for (i = 0; i < nstate; i++) {//Parse all the eigenstates.
			parse_eigenstate(eigls->eigtab[i], eigls->med, i, narg, args);
		}
	}
	sscanf(get_value(narg, args, "realtime"), "%lg", &(eigls->realtime));
}
