/****
 * @date Created on 2021-07-10 at 12:35:29 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to manipulate differential cross section data.
 ***/
#include <stdlib.h>                    /* Standard Library for Memory Allocation */
#include <stdio.h>                     /* Standard Library for Input and Output */
#include <string.h>                    /* Standard Library for String Manipulation */
#include <math.h>                      /* Standard Library of Mathematical Function */
#include "common_util.h"               /* Import the Common Utilities */
#include "green_bessel.h"              /* Import the Bessel Functions */
#include "real_vector_util.h"          /* Import the Real Vector Utilities */
#include "medium_util.h"               /* Import the Random Medium Utilities */
#include "diff_cross_section_type.h"   /* Import the Differential Cross Section Type */
/**
 * Defines some macros:
 */
#define PI          3.1415926535897932385  //Pi constant.
#define SQRTTWOPI   2.5066282746310005024  //Square root of 2*pi.
/**
 * Initializes the random function structure (assuming it can be dereferenced), allocating the data array.
 * The content of the random function must be freed using the corresponding del_*() function.
 */
void init_diff_cross_section(DiffCrossSection* dcs, int nseed, uint64_t iseed, double k, int nbin, int nq, int nc, double thmin, double thmax, char* title) {
	if (nseed <= 0) {
		fprintf(stderr, "[ERROR] Invalid number of seeds '%d', aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	if (nbin <= 2 || nbin > 1e5) {
		fprintf(stderr, "[ERROR] Invalid number of bins '%d', aborting...\n", nbin);
		exit(EXIT_FAILURE);
	}
	if (nq <= 1) {
		fprintf(stderr, "[ERROR] Invalid number of quantiles '%d', aborting...\n", nq);
		exit(EXIT_FAILURE);
	}
	if (nc <= 0 || nc > nseed) {
		fprintf(stderr, "[ERROR] Invalid number of curves '%d', should be in 1:%d, aborting...\n", nc, nseed);
		exit(EXIT_FAILURE);
	}
	if (thmin < 1e-6 || thmin > 180 || thmax < 1e-6 || thmax > 180 || thmin > thmax) {
		fprintf(stderr, "[ERROR] Invalid scattering angles thmin=%g, thmax=%g (expected 1e-6:180), aborting...\n", thmin, thmax);
		exit(EXIT_FAILURE);
	}
	dcs->nseed = nseed;
	dcs->iseed = iseed;
	dcs->k = k;
	dcs->nbin = nbin;
	dcs->nq = nq;
	dcs->ntot = nbin*nq;
	dcs->thmin = thmin;
	dcs->thmax = thmax;
	strncpy(dcs->title, title, 79); //Maximum size of title.
	dcs->data = (double*)calloc(dcs->ntot, sizeof(double));
	dcs->mean = (double*)calloc(dcs->nbin, sizeof(double));
	dcs->nc = nc;
	dcs->curve = (double*)calloc(dcs->nc*dcs->nbin, sizeof(double));
}

/**
 * Deletes the content of the random function, not the pointer itself.
 */
void del_diff_cross_section(DiffCrossSection* dcs) {
	free(dcs->data);
	free(dcs->mean);
	free(dcs->curve);
}

/**
 * Prints a short summary of the parameters used in the random function to the standard output.
 */
void print_param_diff_cross_section(DiffCrossSection* dcs) {
	printf("[INFO] Differential cross section entitled '%s' at k=%gsp^-1 with %d seeds in the angle interval %g:%g using %d bins and %d quantiles.\n",
		dcs->title, dcs->k, dcs->nseed, dcs->thmin, dcs->thmax, dcs->nbin, dcs->nq);
}

/**
 * Returns the scattering angle "theta" (in degrees) corresponding to the bin of given index "i" in the range from i=0 to i=nbin-1.
 */
double get_arg_dcs(DiffCrossSection* dcs, int i) {
	return dcs->thmin*exp(i*log(dcs->thmax/dcs->thmin)/(dcs->nbin-1));
}

/**
 * Returns the Airy differential cross section in the far-field regime for small deflections.
 * The radius of the opaque disk is "a" (units of "sp"), the dimension is "d", and the wave number is "k" (units of sp^-1).
 * The scattering angle "theta" is given in radians from theta=0 to theta=pi.
 */
double airy_pattern_dcs(Medium* med, double k, double theta) {//Private function
	int d = med->d;
	double a = radius_of_ball(med); //Radius of the medium (for ball-shaped gas only).
	double j = 4*creal(free_green_imag(d+1, k*a*SQRTTWOPI, theta/SQRTTWOPI));
	return j*j/pow(k, d-1);
}

/**
 * Returns the differential cross section under the ballistic, i.e., perturbatice, approximation valid for large mean free path (n*sigma*R << 1).
 * The radius of the medium is "a" (units of "sp"), the dimension is "d", and the wave number is "k" (units of sp^-1).
 * The scattering angle "theta" is given in radians from theta=0 to theta=pi.
 */
double ballistic_approx_dcs(Medium* med, double k, double theta) {//Private function
	int n = med->n, d = med->d;
	double a = radius_of_ball(med); //Radius of the medium (for ball-shaped gas only).
	double q = 2*k*sin(theta/2); //Transferred momentum.
	double j = creal(free_green_imag(d+2, q, a)/free_green_imag(d+2, q, 0.)); //Normalized Bessel function.
	double ds = creal(cross_section(med, k))/uball_surface(d); //Differential cross section of a single atom.
	return n*ds*(1. + (n-1)*j*j);
}

/**
 * Computes the minimum and maximum differential cross sections within the quantiles Q(1) and Q(nq-2).
 * This function is only used for plotting purposes.
 */
void find_y_range_dcs(DiffCrossSection* dcs, double* ymin, double* ymax) {//Private function
	*ymin = dcs->data[1];
	*ymax = dcs->data[dcs->nq-2];
	double curymin, curymax;
	int i;
	for (i = 1; i < dcs->nbin; i++) {//Loop on the data.
		curymin = dcs->data[i*dcs->nq + 1];
		curymax = dcs->data[i*dcs->nq + dcs->nq-2];
		if (curymin < *ymin)
			*ymin = curymin;
		if (curymax > *ymax)
			*ymax = curymax;
	}
	//printf("[INFO] Vertical plot bounds: ymin=%g, ymax=%g \n", *ymin, *ymax);
}

/**
 * Plots the Airy diffraction cross section in the TikZ file.
 * This model is valid for large "k" and large number of atoms).
 */
void plot_airy_pattern_dcs(DiffCrossSection* dcs, Medium* med, FILE* tikzfp) {//Private function
	int i;
	double th;
	fprintf(tikzfp, "\\addplot[black, thick, densely dashdotted] coordinates {%%Airy diffraction pattern (%dD)\n\t", med->d-1);
	for (i = 0; i < dcs->nbin; i++) {//Loop on the samples of the Airy pattern diff cross section.
		th = get_arg_dcs(dcs, i); //Abscissa angles in degrees.
		fprintf(tikzfp, "(%.16g, %.16g) ", th, airy_pattern_dcs(med, dcs->k, th*PI/180));
	}
	fprintf(tikzfp, "}; \\addlegendentry{Airy pattern}\n");
}

/**
 * Plots the Airy diffraction cross section in the TikZ file.
 * This model is valid for large "k" and large number of atoms).
 */
void plot_ballistic_dcs(DiffCrossSection* dcs, Medium* med, FILE* tikzfp) {//Private function
	int i;
	double th;
	fprintf(tikzfp, "\\addplot[black, thick, densely dotted] coordinates {%%Ballistic approx\n\t");
	for (i = 0; i < dcs->nbin; i++) {//Loop on the samples of the ballistic diff cross section.
		th = get_arg_dcs(dcs, i); //Abscissa angles in degrees.
		fprintf(tikzfp, "(%.16g, %.16g) ", th, ballistic_approx_dcs(med, dcs->k, th*PI/180));
	}
	fprintf(tikzfp, "}; \\addlegendentry{Ballistic}\n");
}

/**
 * Plots the sample curves for the differential cross section in the given TikZ file.
 */
void plot_curves_dcs(DiffCrossSection* dcs, FILE* tikzfp) {//Private function
	int i, c;
	double th, sat;
	for (c = 0; c < dcs->nc; c++) {//Loop on the curves to plot.
		sat = round(100.*(1 - (double)c/dcs->nc));  //Saturation of the color of the current curve.
		fprintf(tikzfp, "\\addplot[blue!%g, thick, smooth] coordinates {%%Sample curve seed=%lu\n\t", sat, dcs->iseed + c);
		for (i = 0; i < dcs->nbin; i++) {//Loop on the samples of the ballistic diff cross section.
			th = get_arg_dcs(dcs, i); //Abscissa angles in degrees.
			fprintf(tikzfp, "(%.16g, %.16g) ", th, dcs->curve[c + i*dcs->nc]);
		}
		fprintf(tikzfp, "}; \\addlegendentry{$i_{\\rm seed}=%lu$}\n", dcs->iseed + c);
	}
}

/**
 * Plots the main data of the differential cross section in the given TikZ file.
 */
void plot_data_dcs(DiffCrossSection* dcs, FILE* tikzfp) {//Private function
	int i, q, ismedian; //i=Current bin index, q=Current quantile index.
	double th; //Current scattering angle "theta".
	for (q = 1; q < dcs->nq-1; q++) {//Write the data of the quantiles (excluding extreme quantiles q=0 and q=nq-1).
		ismedian = (2*q == dcs->nq-1); //Check if the current quantile is the median.
		if (ismedian) //If the quantile corresponds to the median (only for odd values of "nq").
			fprintf(tikzfp, "\\addplot[blue, thick, densely dashed] coordinates {%%Median curve (quantile %d/%d)\n", q, dcs->nq-1);
		else 
			fprintf(tikzfp, "\\addplot[name path=q%d, draw=none, forget plot] coordinates {%%Quantile %d/%d \n\t", q, q, dcs->nq-1);
		for (i = 0; i < dcs->nbin; i++) {//Loop on the bins.
			th = get_arg_dcs(dcs, i); //Abscissa angles in degrees.
			fprintf(tikzfp, "(%.16g, %.16g) ", th, dcs->data[i*dcs->nq+q]);
		}
		fprintf(tikzfp, "}; ");
		if (ismedian) //Appends a legend entry for the median curve.
			fprintf(tikzfp, "\\addlegendentry{Median curve}");
		fprintf(tikzfp, "\n");
	}
	int cval, qc, qmax = dcs->nq/2; //Maximum quantile for which the interval should be drawn.
	for (q = 1; q < qmax; q++) {//Plot the quantiles (excluding extreme quantiles q=0).
		qc = dcs->nq-q-1; //Complement quantile index.
		cval = (int)round(40.*(1.+q)/qmax); //Color value intensity (0 is white, 100 is pure color).
		fprintf(tikzfp, "\\addplot[blue!%d] fill between[of=q%d and q%d]; \\addlegendentry{$Q_{%d/%d}-Q_{%d/%d}$}\n",
			cval, q, qc, qc, dcs->nq-1, q, dcs->nq-1);
	}
	fprintf(tikzfp, "\\addplot[blue, thick] coordinates {%%Mean curve\n\t"); //Plot the mean curve.
	for (i = 0; i < dcs->nbin; i++) {//Loop on the bins.
		th = get_arg_dcs(dcs, i); //Abscissa angles in degrees.
		fprintf(tikzfp, "(%.16g, %.16g) ", th, dcs->mean[i]);
	}
	fprintf(tikzfp, "}; \\addlegendentry{Mean curve}\n");
}

/**
 * Exports the given differential cross section structure "dcs" in a tikzpicture environment to the file "fname" (given without extension).
 * If "thlog" is 1, then use logarithmic axes for the abscissa, i.e., the scattering angles.
 */
void export_tikz_diff_cross_section(DiffCrossSection* dcs, Medium* med, int thlog, const char* fname) {
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a");
	double ymin, ymax;
	find_y_range_dcs(dcs, &ymin, &ymax); //Determine the vertical range ymin:ymax.
	fprintf(tikzfp, "\\usepgfplotslibrary{fillbetween}\n\\par\\begin{tikzpicture}[%%\n");
	fprintf(tikzfp, "\t/pgfplots/every axis/.style={%%\n\t\twidth=250pt, height=180pt,\n\t\t");
	if (thlog) {
		fprintf(tikzfp, "xmode=log, ");
	}
	fprintf(tikzfp, "ymode=log,\n\t\txmin=%g, xmax=%g,\n\t\tymin=%g, ymax=%g,\n", dcs->thmin, dcs->thmax, ymin/2, 2*ymax);
	fprintf(tikzfp, "\t\txlabel={$\\theta$ ($^\\circ$)},\n\t\tylabel={$\\mathrm{d}\\sigma/\\mathrm{d}\\Omega$ ($\\varsigma^{%d}$)},\n", med->d-1);
	if (thlog) {//If thlog=1, then use log scale for the scattering angle.
		fprintf(tikzfp, "\t\txtick={0.1,0.2,0.5,1,2,5,10,20,45,90,180},\n\t\tminor xtick={0.1,0.2,...,2,3,4,...,20,25,30,...,90,100,110,...,180},\n");
		fprintf(tikzfp, "\t\txticklabel={\\pgfmathparse{exp(\\tick)}\\pgfmathprintnumber[fixed relative, precision=3]{\\pgfmathresult}},\n");
	}
	fprintf(tikzfp, "\t\tyticklabel style={rotate=90},\n\t\ttitle style={align=center},\n\t\tlegend pos=outer north east,\n\t\taxis on top,\n\t},\n]%%\n");
	fprintf(tikzfp, "\\begin{axis}[%%\n\tname={samples},\n");
	fprintf(tikzfp, "\ttitle={%dD %s for $N=%d$ and %s\\\\ at $k=%g\\,\\varsigma^{-1}$ with %d seeds (%.4g s)},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, dcs->k, dcs->nseed, dcs->realtime);
	fprintf(tikzfp, "\txlabel=\\empty,\n\txticklabel=\\empty,\n]%%\n");
	plot_curves_dcs(dcs, tikzfp);
	plot_ballistic_dcs(dcs, med, tikzfp);
	plot_airy_pattern_dcs(dcs, med, tikzfp);
	fprintf(tikzfp, "\\end{axis}%%\n\\begin{axis}[%%\n\tname={stat},\n\tat={(samples.below south west)},\n\tanchor={above north west},\n\tyshift=-3pt,\n]%%\n");
	plot_data_dcs(dcs, tikzfp);
	plot_ballistic_dcs(dcs, med, tikzfp);
	plot_airy_pattern_dcs(dcs, med, tikzfp);
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	fclose(tikzfp);
}

/**
 * Saves the differential cross section data, in addition to other metadata, to the file "fname".
 */
void save_diff_cross_section(DiffCrossSection* dcs, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[diff_cross_section]\ntitle=%s\n", dcs->title);
	fprintf(fp, "nseed=%d\n", dcs->nseed);
	fprintf(fp, "iseed=%lu\n", dcs->iseed);
	fprintf(fp, "k=%g\n", dcs->k);
	fprintf(fp, "thmin=%.16g\n", dcs->thmin);
	fprintf(fp, "thmax=%.16g\n", dcs->thmax);
	fprintf(fp, "nbin=%d\n", dcs->nbin);
	fprintf(fp, "nq=%d\n", dcs->nq);
	fprintf(fp, "realtime=%.16g\n", dcs->realtime);
	fprintf(fp, "data=");
	int l;
	for (l = 0; l < dcs->ntot; l++) {
		fprintf(fp, "%.16g ", dcs->data[l]);
	}
	fprintf(fp, "\nmean=");
	for (l = 0; l < dcs->nbin; l++) {
		fprintf(fp, "%.16g ", dcs->mean[l]);
	}
	fprintf(fp, "\nnc=%d\ncurve=", dcs->nc);
	for (l = 0; l < dcs->nbin*dcs->nc; l++) {
		fprintf(fp, "%.16g ", dcs->curve[l]);
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}

/**
 * Parses the differential cross section structure "dcs" with the arguments of "args".
 * The "dcs" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_diff_cross_section(DiffCrossSection* dcs, int narg, char** args) {
	int nargmin = 8; //Minimum number of argument.
	if (narg < nargmin) {
		fprintf(stderr, "[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	char *value, *title;
	int nbin, nseed, nq, nc;
	uint64_t iseed;
	double k, thmin, thmax;
	value = get_value(narg, args, "nseed");
	if (sscanf(value, "%d", &nseed) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of seeds '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "iseed");
	if (sscanf(value, "%lu", &iseed) != 1) {
		fprintf(stderr, "[ERROR] Invalid initial seed '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "k");
	if (sscanf(value, "%lg", &k) != 1) {
		fprintf(stderr, "[ERROR] Invalid wave number '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "thmin");
	if (sscanf(value, "%lg", &thmin) != 1) {
		fprintf(stderr, "[ERROR] Invalid minimum angle '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "thmax");
	if (sscanf(value, "%lg", &thmax) != 1) {
		fprintf(stderr, "[ERROR] Invalid maximum angle '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nbin");
	if (sscanf(value, "%d", &nbin) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nq");
	if (sscanf(value, "%d", &nq) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of quantiles '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nc");
	if (sscanf(value, "%d", &nc) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of sample curves '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	title = get_value(narg, args, "title");
	init_diff_cross_section(dcs, nseed, iseed, k, nbin, nq, nc, thmin, thmax, title);
	sscanf(get_value(narg, args, "realtime"), "%lg", &(dcs->realtime));
	value = get_value(narg, args, "data"); //Now parses the "data" field.
	if (value[0]) {//If existing "data" field, then parses it.
		int ncount = count_real_data(value);
		if (ncount != dcs->ntot) {
			fprintf(stderr, "[ERROR] Differential cross section data corruption. Expected %d data points but found %d, aborting...\n", dcs->ntot, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(dcs->ntot, dcs->data, value);
	}
	value = get_value(narg, args, "mean"); //Now parses the "mean" field.
	if (value[0]) {//If existing "mean" field, then parses it.
		int ncount = count_real_data(value);
		if (ncount != dcs->nbin) {
			fprintf(stderr, "[ERROR] Differential cross section data corruption. Expected %d mean points but found %d, aborting...\n", dcs->nbin, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(dcs->nbin, dcs->mean, value);
	}
	value = get_value(narg, args, "curve"); //Now parses the "curve" field.
	if (value[0]) {//If existing "curve" field, then parses it.
		int ncount = count_real_data(value);
		if (ncount != dcs->nbin*dcs->nc) {
			fprintf(stderr, "[ERROR] Differential cross section data corruption. Expected %d curve points but found %d, aborting...\n", dcs->nbin*dcs->nc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(dcs->nbin*dcs->nc, dcs->curve, value);
	}
}
