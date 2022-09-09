/****
 * @date Created on 2021-12-21 at 19:25:38 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to create and manipulate longitudinal cuts of the wave function.
 ***/
#include <stdlib.h>                   /* Standard Library for Memory Allocation */
#include <stdio.h>                    /* Standard Library for Input and Output */
#include <string.h>                   /* Standard Library for String Manipulation */
#include <math.h>                     /* Standard Library for Mathematical Functions */
#include "common_util.h"              /* Import the General Purpose Functions */
#include "wave_function_cut_type.h"   /* Import the Wave Function Cut Type */
#include "incident_wave_util.h"       /* Import the Incident Wave Utilities */
#include "real_vector_util.h"         /* Import the Complex Data Parsing functions */
#include "green_bessel.h"             /* Import the Green Functions */
#include "scattering_model_util.h"    /* Import the Scattering Model Utilities */
#include "medium_util.h"              /* Import the Medium Management Utilities */

/**
 * Initialization function of an empty complex map with the given parameters.
 * The content of the complex map must be deleted with the corresponding del_*() function.
 */
void init_wfcut(WaveFunctionCut* wfcut, int nbin, double xmin, double xmax, int nseed, uint64_t iseed, IncidentWave* iwave) {
	if (nbin > 1e7) {//Beyond 10^7 is a bit too large...
		printf("[ERROR] The number of bins is too large (nbin=%d), aborting...\n", nbin);
		exit(EXIT_FAILURE);
	}
	if (nbin < 2) {
		printf("[ERROR] The number of bins is too small (nbin=%d), aborting...\n", nbin);
		exit(EXIT_FAILURE);
	}
	if (xmin >= xmax) {
		printf("[ERROR] Invalid x range (xmin=%g, xmax=%g), aborting...\n", xmin, xmax);
		exit(EXIT_FAILURE);
	}
	if (nseed <= 0) {
		printf("[ERROR] Invalid number of seeds (nseed=%d), aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	wfcut->nbin = nbin;
	wfcut->xmin = xmin;
	wfcut->xmax = xmax;
	wfcut->dx = (xmax - xmin)/(nbin - 1);
	wfcut->iwave = (IncidentWave*)calloc(1, sizeof(IncidentWave));
	copy_iwave(iwave, wfcut->iwave);
	wfcut->nseed = nseed;
	wfcut->iseed = iseed;
	wfcut->realtime = 0.; //Computation time is zero, since empty.
	wfcut->ntot = 3*nbin; //Size of the arrays.
	wfcut->repsi = (double*)calloc(wfcut->ntot, sizeof(double));
	wfcut->impsi = (double*)calloc(wfcut->ntot, sizeof(double));
	wfcut->density = (double*)calloc(wfcut->ntot, sizeof(double));
}

/**
 * Deletes the content (allocated pointers) within the given wave function cut, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_wfcut(WaveFunctionCut* wfcut) {
	free(wfcut->iwave);
	free(wfcut->repsi);
	free(wfcut->impsi);
	free(wfcut->density);
}

/**
 * Prints a short summary of the parameters used in the wave function cut for the standard output.
 */
void print_param_wfcut(WaveFunctionCut* wfcut) {
	printf("[INFO] Wave function with %d bins in region x=%g:%g for %d seeds with source %s.\n",
		wfcut->nbin, wfcut->xmin, wfcut->xmax, wfcut->nseed, wfcut->iwave->desc_ascii);
}

/**
 * Returns the 1D position corresponding to the given index "l" in the range {0, ..., nbin-1}.
 * Note that the other directions in "pos" are unchanged, and typically are left at zero.
 * The samples are assumed to be stored in row-major format in wfcut->function and wfcut->density.
 */
void get_pos_wfcut(WaveFunctionCut* wfcut, int l, double* pos) {
	pos[0] = (wfcut->xmin + l*wfcut->dx);
}

/**
 * Store the given samples data in "wfcut" by computing the average, the first and third quartiles of the samples.
 * This function assumes the samples array are all of size "nseed*nbin", and stored in standard order:
 * samples = {bin1:[seed1, seed2, ...], bin2:[seed1, seed2, ...], bin3:[seed1, seed2, ...], ...}.
 */
void store_data_wfcut(WaveFunctionCut* wfcut, double* repsi_samples, double* impsi_samples, double* density_samples) {
	int i, ns = wfcut->nseed;
	double q1 = 0.25, q3 = 0.75; //Values of the quantiles. There are just 1/4 and 3/4 for the first and third quartiles.
	for (i = 0; i < wfcut->nbin; i++) {//Loop on the bins.
		//Real part of the wave function:
		sort_vector(ns, repsi_samples + i*ns); //The samples must be sorted before computing quantiles (this may take some time).
		wfcut->repsi[3*i]   = mean_vector(ns, repsi_samples + i*ns); //Average value.
		wfcut->repsi[3*i+1] = quantile(ns, repsi_samples + i*ns, q1); //First quartile.
		wfcut->repsi[3*i+2] = quantile(ns, repsi_samples + i*ns, q3); //Third quartile.
		
		//Real part of the wave function:
		sort_vector(ns, impsi_samples + i*ns);
		wfcut->impsi[3*i]   = mean_vector(ns, impsi_samples + i*ns); //Average value.
		wfcut->impsi[3*i+1] = quantile(ns, impsi_samples + i*ns, q1); //First quartile.
		wfcut->impsi[3*i+2] = quantile(ns, impsi_samples + i*ns, q3); //Third quartile.
		
		//Suare modulus of the wave function:
		sort_vector(ns, density_samples + i*ns);
		wfcut->density[3*i]   = mean_vector(ns, density_samples + i*ns); //Average value.
		wfcut->density[3*i+1] = quantile(ns, density_samples + i*ns, q1); //First quartile.
		wfcut->density[3*i+2] = quantile(ns, density_samples + i*ns, q3); //Third quartile.
	}
}

/**
 * Computes and returns the average density for a spherically symmetrical source at radius "r" from the source.
 * "reff" is the value of the effective radius, which is calculated and returned.
 * "amp" is the amplitude of the Dirac delta, which is calculated and returned.
 */
double spherical_density_diffusion(WaveFunctionCut* wfcut, Medium* med, double r, double* reff, double* amp) {//Private
	int d = med->d;
	dcomplex k = wfcut->iwave->k;  //Wavenumber.
	dcomplex i0 = free_green_imag(d, k, 0.); //Value of I(k,0).
	double csec = creal(cross_section(med, k)); //Total cross section of an individual scatterer.
	double surfd = uball_surface(d); //Surface of the unit d-ball.
	double radius = radius_of_ball(med); //Compute the radius of the medium. Only valid in the ball case.
	double delta = 2*uball_volume(d-1)/(surfd*csec); //Variation of the effective radius with respect to the actual radius.
	*amp = d*csec*i0/cabs(k);  //Amplitude of the Dirac delta.
	double density;
	if (d != 2) {//In any dimension different from d=2.
		double p = 2-d; //Power law.
		*reff = radius*pow(1. + p*delta/radius, 1./p);
		density = -(*amp/surfd)*(pow(r, p) - pow(*reff, p))/p;
	}
	else {//Special case d=2.
		*reff = radius*exp(delta/radius);
		density = -(*amp/surfd)*log(r/(*reff));
	}
	return density;
}

/**
 * Appends a TikZ picture environment to a file "fname".
 * Add the metadata in "wfcut" to the plot.
 * Note that the filename "fname" should not have any extension.
 */
void export_tikz_wfcut(WaveFunctionCut* wfcut, Medium* med, const char* fname) {//Public
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	int l, d = med->d;
	double* pos = (double*)calloc(d, sizeof(double)); //Current position.
	double radius = radius_of_ball(med); //Compute the radius of the medium. Only valid in the ball case.
	//\t/pgfplots/extra x ticks={%.16lg,%.16lg},\n\t/pgfplots/extra x tick labels={$-R$,$R$},\n\t/pgfplots/extra x tick style={grid=major, tick label style={rotate=90, anchor=south west}},\n
	fprintf(tikzfp, "\\usepgfplotslibrary{fillbetween}%%\n\\par%%\n\\begin{tikzpicture}[%%\n\t/pgfplots/width=250pt,\n\t/pgfplots/xlabel={$x~(\\varsigma)$},\n\t/pgfplots/title style={align=center},\n\t/pgfplots/legend pos={outer north east},\n\t/pgfplots/enlargelimits=false,\n\t/pgfplots/yticklabel style={rotate=90},\n\t/pgfplots/axis on top=true,\n]%%\n\\pgfplotstableread{%% Raw wave function data.\n"); //, -radius, +radius
	dcomplex gavg, gfree; //Local values of Green functions.
	dcomplex k = wfcut->iwave->k, keff = csqrt(k*k - 1./invf_scmodel(med->scmodel, d, k)); //Effective wavenumber keff=sqrt(k^2-nF(k)).
	fprintf(tikzfp, "\tx\trepsi_avg\trepsi_q1\trepsi_q3\tgreen_avg\tgreen_free\n");
	for (l = 0; l < wfcut->nbin; l++) {//Loop on the samples of the mean curve.
		get_pos_wfcut(wfcut, l, pos);
		gavg = free_green(d, keff, norm(d, pos));
		gfree = free_green(d, k, norm(d, pos));
		fprintf(tikzfp, "\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
			pos[0], wfcut->repsi[3*l], wfcut->repsi[3*l+1], wfcut->repsi[3*l+2], creal(gavg), creal(gfree));
	}
	fprintf(tikzfp, "}\\repsidata%%\n\\begin{axis}[%%\n\tname={function},\n\ttitle={%dD %s for $N=%d$ and %s\\\\ $\\langle\\psi(\\mathbf{r})\\rangle$, source %s\\\\ with %d seeds (%d bins, %.1f s)},\n\tylabel={$\\psi(\\mathbf{r})$},\n]%%\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, wfcut->iwave->desc_latex, wfcut->nseed, wfcut->nbin, wfcut->realtime);
	fprintf(tikzfp, "\\addplot[blue, thick] table[y=repsi_avg] {\\repsidata}; \\addlegendentry{$\\langle\\Re\\psi(\\mathbf{r})\\rangle$}\n");
	fprintf(tikzfp, "\\addplot[name path=req1, draw=none, forget plot] table[y=repsi_q1] {\\repsidata};\n");
	fprintf(tikzfp, "\\addplot[name path=req3, draw=none, forget plot] table[y=repsi_q3] {\\repsidata};\n");
	fprintf(tikzfp, "\\addplot[blue, opacity=0.4] fill between[of=req1 and req3]; \\addlegendentry{Interquartile}\n");
	//fprintf(tikzfp, "\\addplot[red, thick] table[y=impsi_avg] {\\psidata}; \\addlegendentry{$\\langle\\Im\\psi(\\mathbf{r})\\rangle$}\n");
	//fprintf(tikzfp, "\\addplot[name path=imq1, draw=none, forget plot] table[y=impsi_q1] {\\psidata};\n");
	//fprintf(tikzfp, "\\addplot[name path=imq3, draw=none, forget plot] table[y=impsi_q3] {\\psidata};\n");
	//fprintf(tikzfp, "\\addplot[red, opacity=0.4] fill between[of=imq1 and imq3]; \\addlegendentry{Interquartile}\n");
	if (wfcut->iwave->source == spherical) {//In case of spherical source, then show the theoretical average Green function.
		fprintf(tikzfp, "\\addplot[black, thick, densely dashed] table[y=green_avg] {\\repsidata}; \\addlegendentry{$\\Re G^+(\\kappa,r)$} %% Avg Green fct, keff=%g%+gi\n", creal(keff), cimag(keff));
		fprintf(tikzfp, "\\addplot[black, thick, densely dotted] table[y=green_free] {\\repsidata}; \\addlegendentry{$\\Re G^+(k,r)$} %% Free Green fct, k=%g%+gi\n", creal(k), cimag(k));
	}
	fprintf(tikzfp, "\\draw[black!60] (axis cs:%.15lg,\\pgfkeysvalueof{/pgfplots/ymin}) -- (axis cs:%.15lg,\\pgfkeysvalueof{/pgfplots/ymax}) node[pos=0,above right,rotate=90]{$R$}; %%Radius of the medium.\n", radius, radius);
	fprintf(tikzfp, "\\end{axis}%%\n");
	fprintf(tikzfp, "\\pgfplotstableread{%% Raw density data.\n");
	fprintf(tikzfp, "\tx\tdensity_avg\tdensity_q1\tdensity_q3\tdiffusion\tgreen_avg_abs2\n");
	double density_diff, reff, amp;
	for (l = 0; l < wfcut->nbin; l++) {//Loop on the samples of the mean curve.
		get_pos_wfcut(wfcut, l, pos);
		density_diff = spherical_density_diffusion(wfcut, med, norm(d, pos), &reff, &amp); //Theoretical density in the diffusion approximation.
		gavg = free_green(d, keff, norm(d, pos));
		fprintf(tikzfp, "\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\t%.16lg\n",
			pos[0], wfcut->density[3*l], wfcut->density[3*l+1], wfcut->density[3*l+2], density_diff, creal(gavg)*creal(gavg) + cimag(gavg)*cimag(gavg));
	}
	fprintf(tikzfp, "}\\densitydata%%\n\\begin{axis}[%%\n\tname={density},\n\tat={(function.below south west)},\n\tanchor={above north west},\n\tyshift=-10pt,\n\tymode=log,\n\ttitle={%dD %s for $N=%d$ and %s\\\\ $\\langle|\\psi(\\mathbf{r})|^2\\rangle$, source %s\\\\ with %d seeds (%d bins, %.1f s)},\n]%%\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, wfcut->iwave->desc_latex, wfcut->nseed, wfcut->nbin, wfcut->realtime);
	fprintf(tikzfp, "\\addplot[blue, thick] table[y=density_avg] {\\densitydata}; \\addlegendentry{$\\langle|\\psi(\\mathbf{r})|^2\\rangle$}\n");
	fprintf(tikzfp, "\\addplot[name path=absq1, draw=none, forget plot] table[y=density_q1] {\\densitydata};\n");
	fprintf(tikzfp, "\\addplot[name path=absq3, draw=none, forget plot] table[y=density_q3] {\\densitydata};\n");
	fprintf(tikzfp, "\\addplot[blue, opacity=0.4] fill between[of=absq1 and absq3]; \\addlegendentry{Interquartile}\n");
	if (wfcut->iwave->source == spherical) {//In case of spherical source, then show the theoretical average Green function.
		fprintf(tikzfp, "\\addplot[black, thick, densely dashed] table[y=diffusion] {\\densitydata}; \\addlegendentry{$\\rho_{\\rm da}(r)$} %% Density in diffusion approx (da), reff=%g, amp=%g\n", reff, amp);
		fprintf(tikzfp, "\\addplot[black, thick, densely dotted] table[y=green_avg_abs2] {\\densitydata}; \\addlegendentry{$|G^+(\\kappa(k),r)|^2$} %% Avg Green fct, keff=%g%+gi\n", creal(keff), cimag(keff));
	}
	fprintf(tikzfp, "\\draw[black!60] (axis cs:%.15lg,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%.15lg,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[pos=0,above right,rotate=90]{$R$}; %%Radius of the medium.\n", radius, radius);
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	fclose(tikzfp);
}

/**
 * Saves the wave function cut data, in addition to other metadata, to the file titled "fname".
 */
void save_wfcut(WaveFunctionCut* wfcut, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[wavefunction_cut]\n");
	fprintf(fp, "nseed=%d\n", wfcut->nseed);
	fprintf(fp, "iseed=%lu\n", wfcut->iseed);
	fprintf(fp, "nbin=%d\n", wfcut->nbin);
	fprintf(fp, "xrange=%.16lg:%.16lg\n", wfcut->xmin, wfcut->xmax);
	save_iwave(wfcut->iwave, fp);
	fprintf(fp, "realtime=%g\n", wfcut->realtime); //Total physical duration of the computation of the complex map (in seconds).
	save_real_data(wfcut->ntot, wfcut->repsi, "repsi", fp);
	save_real_data(wfcut->ntot, wfcut->impsi, "impsi", fp);
	save_real_data(wfcut->ntot, wfcut->density, "density", fp);
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Parses the wave function cut "wfcut" with the arguments of "args".
 * The "wfcut" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_wfcut(WaveFunctionCut* wfcut, int narg, char** args) {
	int nargmin = 6; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	int nbin = 0;  //Default value.
	int nseed;
	double xmin, xmax;
	uint64_t iseed;
	char *value;
	value = get_value(narg, args, "nbin");
	if (sscanf(value, "%d", &nbin) != 1) {
		printf("[ERROR] Invalid number of sample nbin='%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "xrange");
	if (sscanf(value, "%lg:%lg", &xmin, &xmax) != 2) {
		printf("[ERROR] Invalid x range '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nseed");
	if (sscanf(value, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "iseed");
	if (sscanf(value, "%lu", &iseed) != 1) {
		printf("[ERROR] Invalid initial seed '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	IncidentWave iwave; //Parse the incident wave.
	if (parse_iwave(&iwave, narg, args)) {//If the incident wave parsing failed
		printf("[ERROR] Invalid incident wave, aborting...\n");
		exit(EXIT_FAILURE);
	}
	init_wfcut(wfcut, nbin, xmin, xmax, nseed, iseed, &iwave);
	sscanf(get_value(narg, args, "realtime"), "%lg", &(wfcut->realtime)); //Set realtime after initialization.
	int ncount, ncount_expc = wfcut->ntot; //Expected number of real points.
	value = get_value(narg, args, "repsi");  //Now parsing the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != ncount_expc) {
			fprintf(stderr, "[ERROR] Data corruption in 'repsi'. Expected %d data points but found %d, aborting...\n", ncount_expc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ncount, wfcut->repsi, value);
	}
	value = get_value(narg, args, "impsi");  //Now parsing the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != ncount_expc) {
			fprintf(stderr, "[ERROR] Data corruption in 'impsi'. Expected %d data points but found %d, aborting...\n", ncount_expc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ncount, wfcut->impsi, value);
	}
	value = get_value(narg, args, "density");  //Now parsing the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != ncount_expc) {
			fprintf(stderr, "[ERROR] Data corruption in 'density'. Expected %d data points but found %d, aborting...\n", ncount_expc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ncount, wfcut->density, value);
	}
}

