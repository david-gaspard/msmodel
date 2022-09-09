/****
 * @date Created on 2021-07-05 at 19:18:18 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to manipulate random function structures.
 ***/
#include <stdlib.h>                 /* Standard Library for Memory Allocation */
#include <stdio.h>                  /* Standard Library for Input and Output */
#include <string.h>                 /* Standard Library for String Manipulation */
#include <math.h>                   /* Standard Library for Mathematical Functions */
#include "common_util.h"            /* Import the Common Utilities */
#include "real_vector_util.h"       /* Import the Real Vector Utilities */
#include "medium_util.h"            /* Import the Random Medium Utilities */
#include "random_function_type.h"   /* Import the Random Function Type */
/**
 * Define some macros
 */
#define PI   3.1415926535897932385
/**
 * Initializes the random function structure (assuming it can be dereferenced), allocating the data array.
 * The content of the random function must be freed using the corresponding del_*() function.
 */
void init_random_function(RandomFunction* rfun, RandomFunctionType type, int nseed, int nbin, double xmin, double xmax, char* title) {
	if (nbin <= 2 || nbin > 1e5) {
		fprintf(stderr, "[ERROR] Invalid number of bins '%d', aborting...\n", nbin);
		exit(EXIT_FAILURE);
	}
	if (nseed <= 0) {
		fprintf(stderr, "[ERROR] Invalid number of seeds '%d', aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	if (xmin < 1e-6 || xmax < 1e-6 || xmin > xmax) {
		fprintf(stderr, "[ERROR] Invalid abscissa xmin=%g, xmax=%g, aborting...\n", xmin, xmax);
		exit(EXIT_FAILURE);
	}
	rfun->type = type;
	rfun->nbin = nbin;
	rfun->nseed = nseed;
	rfun->ntot = nbin*nseed;
	rfun->xmin = xmin;
	rfun->xmax = xmax;
	strncpy(rfun->title, title, 79); //Maximum size of title.
	rfun->data = (double*)calloc(rfun->ntot, sizeof(double));
}

/**
 * Deletes the content of the random function, not the pointer itself.
 */
void del_random_function(RandomFunction* rfun) {
	free(rfun->data);
}

/**
 * Converts the given type of random function to the corresponding string.
 * Returns the string "none" if the shape is not recognized. 
 */
char* random_function_type_to_string(RandomFunctionType type) {
	switch (type) {
		case crsec:
			return "crsec";
		case dos:
			return "dos";
		default:
			return "none";
	}
}

/**
 * Converts the given string to the corresponding type of random function.
 * Returns 1 if the string is not recognized. 
 */
int string_to_random_function_type(const char* str, RandomFunctionType* type) {
	if (strcmp(str, "crsec") == 0) {
		*type = crsec;
	}
	else if (strcmp(str, "dos") == 0) {
		*type = dos;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Prints a short summary of the parameters used in the random function to the standard output.
 */
void print_param_random_function(RandomFunction* rfun) {
	printf("[INFO] Random function entitled '%s' of type '%s' with %d bins and %d seeds in the interval %g:%g.\n",
		rfun->title, random_function_type_to_string(rfun->type), rfun->nbin, rfun->nseed, rfun->xmin, rfun->xmax);
}

/**
 * Returns the abscissa corresponding to the bin of given index "i" in the range from i=0 to i=nbin-1.
 */
double get_arg_rfun(RandomFunction* rfun, int i) {
	//return rfun->xmin + i*(rfun->xmax - rfun->xmin)/(rfun->nbin - 1); //The abscissa is the wave number "k" in units of 1/sp.
	return rfun->xmin*exp(i*log(rfun->xmax/rfun->xmin)/(rfun->nbin-1));
}

/**
 * Sorts the samples of each bin in increasing order so that the quantiles can be computed.
 * This function should be called once before saving the data points to a file.
 */
void sort_random_function(RandomFunction* rfun) {
	int i;
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
		sort_vector(rfun->nseed, rfun->data + i*rfun->nseed);
	}
}

/**
 * Exports the given random function structure "rfun" in a tikzpicture environment to the file "fname" (given without extension).
 */
void export_tikz_random_function(RandomFunction* rfun, Medium* med, const char* fname) {
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a");
	fprintf(tikzfp, "\\usepgfplotslibrary{fillbetween}\n\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\txmode=log,\n\txmin=%g, xmax=%g,\n", rfun->xmin, rfun->xmax); //\tymin=0,\n
	fprintf(tikzfp, "\ttitle={%dD %s for $N=%d$ and %s\\\\ with %d bins and %d seeds (%.4g s)},\n\ttitle style={align=center},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, rfun->nbin, rfun->nseed, rfun->realtime);
	fprintf(tikzfp, "\txlabel={$k$ ($\\varsigma^{-1}$)},\n\tylabel={%s},\n", random_function_type_to_string(rfun->type));
	fprintf(tikzfp, "\tlegend pos={outer north east},\n\tyticklabel style={rotate=90},\n\taxis on top,\n]%%\n");
	//Prepare for the plots with the pre-computation of the quantiles:
	int i, nc = 5; //Number of columns of data.
	if (rfun->type == dos) nc++; //Add one curve for the density of states.
	double* plotdata = (double*)calloc(nc*rfun->nbin, sizeof(double)); //Allocate some place for the pre-computed quantiles.
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins to pre-compute the quantiles.
		plotdata[i*nc] = get_arg_rfun(rfun, i); //The abscissa is the wave number "k" in units of sp^-1.
		plotdata[i*nc + 1] = quantile(rfun->nseed, rfun->data + i*rfun->nseed, 0.25); //NB: Also possible to plot interdecile range (0.10 to 0.90).
		plotdata[i*nc + 2] = quantile(rfun->nseed, rfun->data + i*rfun->nseed, 0.50);
		plotdata[i*nc + 3] = quantile(rfun->nseed, rfun->data + i*rfun->nseed, 0.75);
		plotdata[i*nc + 4] = mean_vector(rfun->nseed, rfun->data + i*rfun->nseed); //Compute the mean curve.
		if (rfun->type == dos) {//Semi-classical approx of independent scatterers. WARN: This perturbative calculation is not expected to work with many atoms.
			plotdata[i*nc + 5] = -med->n*med->scmodel->p[0]/(2*PI*plotdata[i*nc]); //Semi-classical approx in 1D and 3D only.
		}
	}
	fprintf(tikzfp, "\\addplot[name path=q1, draw=none, forget plot] coordinates {%%Quantile 1/4 \n\t");
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
		fprintf(tikzfp, "(%.16g, %.16g) ", plotdata[i*nc], plotdata[i*nc+1]);
	}
	fprintf(tikzfp, "\n};\n\\addplot[name path=q3, draw=none, forget plot] coordinates {%%Quantile 3/4 \n\t");
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
		fprintf(tikzfp, "(%.16g, %.16g) ", plotdata[i*nc], plotdata[i*nc+3]);
	}
	fprintf(tikzfp, "\n};\n\\addplot[blue!40] fill between[of=q1 and q3]; \\addlegendentry{Interquartile}\n");
	fprintf(tikzfp, "\\addplot[blue, thick, densely dashed] coordinates {%%Median curve (quantile 2/4)\n\t");
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
		fprintf(tikzfp, "(%.16g, %.16g) ", plotdata[i*nc], plotdata[i*nc+2]);
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{Median curve}\n");
	fprintf(tikzfp, "\\addplot[blue, thick] coordinates {%%Mean curve\n\t");
	for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
		fprintf(tikzfp, "(%.16g, %.16g) ", plotdata[i*nc], plotdata[i*nc+4]);
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{Mean curve}\n");
	if (rfun->type == dos) {
		fprintf(tikzfp, "\\addplot[black, thick, densely dotted] coordinates {%%Semi-classical approx with independent scatterers\n\t");
		for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
			fprintf(tikzfp, "(%.16g, %.16g) ", plotdata[i*nc], plotdata[i*nc+5]);
		}
		fprintf(tikzfp, "\n}; \\addlegendentry{$-N\\alpha/(2\\pi k)$}\n");
	}
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	free(plotdata);
	fclose(tikzfp);
}

/**
 * Saves the random function data, in addition to other metadata, to the file "fname".
 */
void save_random_function(RandomFunction* rfun, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[random_function]\ntitle=%s\n", rfun->title);
	fprintf(fp, "type=%s\n", random_function_type_to_string(rfun->type));
	fprintf(fp, "nseed=%d\n", rfun->nseed);
	fprintf(fp, "nbin=%d\n", rfun->nbin);
	fprintf(fp, "xmin=%.16g\n", rfun->xmin);
	fprintf(fp, "xmax=%.16g\n", rfun->xmax);
	fprintf(fp, "realtime=%.16g\n", rfun->realtime);
	fprintf(fp, "data=");
	int l;
	for (l = 0; l < rfun->ntot; l++) {
		fprintf(fp, "%.16g ", rfun->data[l]);
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}

/**
 * Parses the random function "rfun" with the arguments of "args".
 * The "rfun" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_random_function(RandomFunction* rfun, int narg, char** args) {
	int nargmin = 5; //Minimum number of argument.
	if (narg < nargmin) {
		fprintf(stderr, "[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	RandomFunctionType type;
	char *value, *title;
	int nbin, nseed;
	double xmin, xmax;
	value = get_value(narg, args, "nseed");
	if (sscanf(value, "%d", &nseed) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of seeds '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nbin");
	if (sscanf(value, "%d", &nbin) != 1) {
		fprintf(stderr, "[ERROR] Invalid number of bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "xmin");
	if (sscanf(value, "%lg", &xmin) != 1) {
		fprintf(stderr, "[ERROR] Invalid minimum abscissa '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "xmax");
	if (sscanf(value, "%lg", &xmax) != 1) {
		fprintf(stderr, "[ERROR] Invalid maximum abscissa '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "type");
	if (string_to_random_function_type(value, &type)) {
		fprintf(stderr, "[ERROR] Invalid random function type '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	title = get_value(narg, args, "title");
	init_random_function(rfun, type, nseed, nbin, xmin, xmax, title);
	sscanf(get_value(narg, args, "realtime"), "%lg", &(rfun->realtime));
	value = get_value(narg, args, "data"); //Now parses the data field.
	if (value[0]) {//If existing data field, then parses it.
		int ncount = count_real_data(value);
		if (ncount != rfun->ntot) {
			fprintf(stderr, "[ERROR] Random function data corruption. Expected %d data points but found %d, aborting...\n", rfun->ntot, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(rfun->ntot, rfun->data, value);
	}
}

