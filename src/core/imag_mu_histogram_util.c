/****
 * @date Created on 2021-07-31 at 14:10:48 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to create and manipulate histograms of the imaginary part of the "mu" eigenvalues.
 * As a reminder, the complex "mu" eigenvalues are computed from the complex multiple scattering matrix.
 * The imaginary part of "mu" is positive and is distributed approximately as the Marchenko-Pastur law for positive-definite matrices.
 * This histogram uses a logarithmic scale.
 ***/
#include <stdlib.h>                  /* Standard Library for Memory Allocation */
#include <stdio.h>                   /* Standard Library for Input and Output */
#include <string.h>                  /* Standard Library for String Manipulation */
#include <math.h>                    /* Standard Library for Mathematical Functions */
#include "dcomplex_type.h"           /* Import the Complex Type and Functions */
#include "common_util.h"             /* Import the General Purpose Functions */
#include "real_vector_util.h"        /* Import the Real Vectors Utilities */
#include "ms_matrix_util.h"          /* Import the Multiple Scattering Matrix Utilities (used to setup the histogram bounds) */
#include "medium_util.h"             /* Import the Medium Utilities */
#include "imag_mu_histogram_type.h"  /* Import the Type of the Histogram of the Imaginary Part of Mu */
/**
 * Define some macros:
 */
#define PI   3.1415926535897932385   //The pi constant.
/**
 * Initialization function of an empty histogram of the imaginary part of "mu" with the given parameters.
 * The content of the histogram must be deleted with the corresponding del_*() function.
 */
void init_ihist(ImagMuHistogram* ihist, int nbin, int nseed, double k, const char* title) {
	if (nbin < 2 || nbin > 1e6) {
		printf("[ERROR] Invalid number of bins (nbin=%d), aborting...\n", nbin);
		exit(EXIT_FAILURE);
	}
	if (nseed < 1) {
		printf("[ERROR] Invalid number of seeds (nseed=%d), aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	if (k <= 0.0) {
		printf("[ERROR] Invalid wave number (k=%g), aborting...\n", k);
		exit(EXIT_FAILURE);
	}
	ihist->nbin = nbin;
	ihist->nseed = nseed;
	ihist->k = k;
	strncpy(ihist->title, title, 79); //Maximum size of title.
	ihist->data = (double*)calloc(ihist->nbin, sizeof(double)); //The data should be initialized to zero everywhere.
}

/**
 * Deletes the content withing the given histogram, not the pointer itself.
 * Do not forget to invoke after use of "init_*" functions.
 */
void del_ihist(ImagMuHistogram* ihist) {
	free(ihist->data);
}

/**
 * Prints a short summary of the parameters used in the histogram for the standard output.
 */
void print_param_ihist(ImagMuHistogram* ihist) {
	printf("[INFO] Histogram of Im(mu) entitled '%s' with %d bins for k=%g and %d seeds.\n", ihist->title, ihist->nbin, ihist->k, ihist->nseed);
}

/**
 * Returns the position "x" of the center of the i-th bin of the given histogram in the interval {0,..., nbin-1}.
 * Note that the histogram uses a logarithmic scale.
 */
double get_arg_ihist(ImagMuHistogram* ihist, int i) {
	return ihist->xmin*exp(log(ihist->xmax/ihist->xmin)*(i + 0.5)/ihist->nbin);
}

/**
 * Initializes the domain of the histogram using the first seed of the medium.
 * This function must be called before appending data to the histogram.
 */
void setup_domain_ihist(ImagMuHistogram* ihist, Medium* med) {
	uint64_t seed0 = 1.;  //Default seed used to initialize the histogram.
	dcomplex* matrix = (dcomplex*)calloc(med->n*med->n, sizeof(dcomplex));  //Multiple scattering matrix.
	dcomplex* mu = (dcomplex*)calloc(med->n, sizeof(dcomplex));  //List of complex eigenvalues of M.
	fill_medium(med, seed0);   //Fill the medium with the given seed.
	build_ms_matrix(med, ihist->k, matrix);  //Builds the matrix with the current value of k.
	eigvals(med->n, matrix, mu);  //Computes the eigenvalues using LAPACK on private slots "mu". Most time-consuming operation.
	double x, xmin = cimag(mu[0]), xmax = cimag(mu[0]); //Minimum and maximum samples.
	int i;
	for (i = 1; i < med->n; i++) {//Find the smallest and largest values of Im(mu).
		x = cimag(mu[i]);
		if (x < xmin && x > 0.) //Find the smallest positive sample, ignoring spurious negative ones.
			xmin = x;
		if (x > xmax)
			xmax = x;
	}
	//printf("[INFO] Setup domain with xmin=%g, xmax=%g \n", xmin, xmax); fflush(stdout);
	double fac = 2.; //Safety margin factor (>1).
	ihist->xmin = xmin/fac;
	ihist->xmax = fac*xmax;
	free(matrix);
	free(mu);
}

/**
 * Appends a list of samples to the histogram. Returns the number of lost samples.
 * Core function of the histogram object.
 */
int append_list_ihist(ImagMuHistogram* ihist, int nz, dcomplex* z) {
	double logx;
	double logxmin = log(ihist->xmin);
	double logxmax = log(ihist->xmax);
	int i, idx, lost = 0;
	for (i = 0; i < nz; i++) {//Loop on the samples to account for.
		logx = log(cimag(z[i])); //Already takes the logarithm here. This assumes positive imaginary part only.
		idx = (int)floor(ihist->nbin*(logx - logxmin)/(logxmax - logxmin));
		if (0 <= idx && idx < ihist->nbin) //Accounts for the sample only in the interval.
			ihist->data[idx] += 1.;
		else //If not in the interval, then it is lost.
			lost++;
	}
	return lost;
}

/**
 * Computes the statistical moment of order "nu" of the given histogram using numerical integration (midpoint method).
 * This function accounts for the logarithmic scale of the histogram using the Jacobian "x".
 */
double moment_ihist(ImagMuHistogram* ihist, double nu) {
	double x, sum = 0.;
	int i;
	for (i = 0; i < ihist->nbin; i++) {//Loop on the bins.
		x = get_arg_ihist(ihist, i); //Coordinate of the center of the i-th bin.
		sum += ihist->data[i]*pow(x, nu)*x;
	}
	return sum*log(ihist->xmax/ihist->xmin)/ihist->nbin;
}

/**
 * Normalize the given histogram to make its integral being equal to 1. The logarithmic bin size is also accounted for.
 * This function must be called once immediately after finishing to append data to the histogram, and before any other operation.
 */
void normalize_ihist(ImagMuHistogram* ihist) {
	double total = 0.; //Total number of samples.
	int i;
	for (i = 0; i < ihist->nbin; i++) {//Loop on the bins.
		total += ihist->data[i];
	}
	double x, fac = total*log(ihist->xmax/ihist->xmin)/ihist->nbin;
	for (i = 0; i < ihist->nbin; i++) {//Normalize each bin.
		x = get_arg_ihist(ihist, i); //Coordinate of the center of the i-th bin.
		ihist->data[i] /= x*fac; //Each bin must be divided by its width.
	}
}

/**
 * Exports the given histogram to a TikZ file of given name "fname" (without extension).
 */
void export_tikz_ihist(ImagMuHistogram* ihist, Medium* med, const char* fname) {
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a");
	fprintf(tikzfp, "\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\txmode=log, ymode=log,\n\txmin=%g, xmax=%g,\n", ihist->xmin, ihist->xmax); //\txmode=log, ymode=log,\n
	fprintf(tikzfp, "\ttitle={%dD %s for $N=%d$ and %s\\\\ at $k=%g\\,\\varsigma^{-1}$ with %d seeds (%.4g s)},\n\ttitle style={align=center},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, ihist->k, ihist->nseed, ihist->realtime);
	fprintf(tikzfp, "\txlabel={$\\Im\\mu$},\n\tylabel={$p(\\Im\\mu)$},\n");
	fprintf(tikzfp, "\tyticklabel style={rotate=90},\n\tlegend pos={outer north east},\n\taxis on top=true,\n]%%\n");
	fprintf(tikzfp, "\\addplot[blue, thick] coordinates {%%Histogram of Im(mu), avgx=%g, varx=%g\n", ihist->avgx, ihist->varx_true);
	double x;
	int i;
	for (i = 0; i < ihist->nbin; i++) {//Loop on the samples to plot.
		x = get_arg_ihist(ihist, i); //Get the center of the ith bin.
		fprintf(tikzfp, "\t(%g, %g)\n", x, ihist->data[i]);
	}
	fprintf(tikzfp, "}; \\addlegendentry{Histogram}\n");
	double sigmax = sqrt(ihist->varx_expc); //Use the theoretical variance, so that the curve is a prediction.
	double xp = (ihist->avgx + sigmax)*(ihist->avgx + sigmax)/ihist->avgx; //Estimates the Marchenko-Pastur parameters from the theoretical mean and variance.
	double xm = (ihist->avgx - sigmax)*(ihist->avgx - sigmax)/ihist->avgx;
	double cn = 0.5*PI*pow(sqrt(xp)-sqrt(xm), 2); //Normalization coefficient located in the denominator.
	fprintf(tikzfp, "\\addplot[black, densely dotted, thick] coordinates {%%Marchenko-Pastur distribution (xm=%g, xp=%g)\n", xm, xp);
	double px, dlx = log(xp/xm)/ihist->nbin;
	for (i = 0; i < ihist->nbin; i++) {//Loop on the samples to plot.
		x = xm*exp((i+0.5)*dlx); //Position of the center of the i-th interval.
		px = sqrt((xp-x)*(x-xm))/(cn*x); //Marchenko-Pastur density itself.
		fprintf(tikzfp, "\t(%g, %g)\n", x, px);
	}
	fprintf(tikzfp, "}; \\addlegendentry{Marchenko-Pastur}\n");
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	fclose(tikzfp);
}

/**
 * Saves the complex histogram data, in addition to other metadata, to the file "fname".
 */
void save_ihist(ImagMuHistogram* ihist, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[imag_mu]\ntitle=%s\n", ihist->title);
	fprintf(fp, "k=%g\n", ihist->k);
	fprintf(fp, "nseed=%d\n", ihist->nseed);
	fprintf(fp, "nbin=%d\n", ihist->nbin);
	fprintf(fp, "xmin=%.16g\n", ihist->xmin);
	fprintf(fp, "xmax=%.16g\n", ihist->xmax);
	fprintf(fp, "avgx=%.16g\n", ihist->avgx);
	fprintf(fp, "varx_true=%.16g\n", ihist->varx_true);
	fprintf(fp, "varx_expc=%.16g\n", ihist->varx_expc);
	fprintf(fp, "realtime=%g\n", ihist->realtime);
	fprintf(fp, "data=");
	int l;
	for (l = 0; l < ihist->nbin; l++) {
		fprintf(fp, "%.16g ", ihist->data[l]);
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}

/**
 * Parses the given histogram with the arguments of "args".
 * The "ihist" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_ihist(ImagMuHistogram* ihist, int narg, char** args) {
	int nargmin = 3; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	char *value, *title;
	int nbin, nseed;
	double k;
	value = get_value(narg, args, "nbin");
	if (sscanf(value, "%d", &nbin) != 1) {
		printf("[ERROR] Invalid number of bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nseed");
	if (sscanf(value, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "k");
	if (sscanf(value, "%lg", &k) != 1) {
		printf("[ERROR] Invalid wave number '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	title = get_value(narg, args, "title");
	init_ihist(ihist, nbin, nseed, k, title); //Call initialization and allocation function.
	sscanf(get_value(narg, args, "xmin"), "%lg", &(ihist->xmin));
	sscanf(get_value(narg, args, "xmax"), "%lg", &(ihist->xmax));
	sscanf(get_value(narg, args, "avgx"), "%lg", &(ihist->avgx));
	sscanf(get_value(narg, args, "varx_true"), "%lg", &(ihist->varx_true));
	sscanf(get_value(narg, args, "varx_expc"), "%lg", &(ihist->varx_expc));
	sscanf(get_value(narg, args, "realtime"), "%lg", &(ihist->realtime));
	value = get_value(narg, args, "data"); //Now parses the data field.
	if (value[0]) {//If existing data field, then parses it.
		int ncount = count_real_data(value);
		if (ncount != ihist->nbin) {
			fprintf(stderr, "[ERROR] Histogram data corruption. Expected %d data points but found %d, aborting...\n", ihist->nbin, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ihist->nbin, ihist->data, value);
	}
}
