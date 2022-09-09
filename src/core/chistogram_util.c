/****
 * @date Created on 2021-09-14 at 15:59:32 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the complex histogram utilites.
 ***/
#include <stdlib.h>             /* Standard Library for Memory Allocation */
#include <stdio.h>              /* Standard Library for Input and Output */
#include <string.h>             /* Standard Library for String Manipulation */
#include <math.h>               /* Standard Library for Mathematical Functions */
#include "common_util.h"        /* Import the General Purpose Functions */
#include "green_bessel.h"       /* Import the Green and Bessel Functions */
#include "real_vector_util.h"   /* Import the Real Vector Manipulation Utilities */
#include "domain_util.h"        /* Import the Domain Utilities */
#include "color_rule_util.h"    /* Import the Color Rule Utilities */
#include "medium_util.h"        /* Import the Medium Utilities */
#include "chistogram_type.h"    /* Import the Complex Histogram Type */
/**
 * Define some macros:
 */
#define PI   3.1415926535897932385   //The pi constant.
/**
 * Initialization function of an empty complex histogram with the given parameters.
 * The content of the complex histogram must be deleted with the corresponding del_*() function.
 */
void init_chistogram(Chistogram* chist, MatrixType mtype, dcomplex k, int nseed, int xbin, int ybin, int mbin, const char* title) {
	if (xbin < 2 || ybin < 2 || mbin < 2) {
		printf("[ERROR] The number of bins is too small (xbin=%d, ybin=%d, mbin=%d), aborting...\n", xbin, ybin, mbin);
		exit(EXIT_FAILURE);
	}
	if (xbin*ybin > 1e6 || mbin > 1e6) {//Beyond 1000x1000 is a bit too large...
		printf("[ERROR] The number of bins is too large (xbin=%d, ybin=%d, mbin=%d), aborting...\n", xbin, ybin, mbin);
		exit(EXIT_FAILURE);
	}
	if (nseed < 1) {
		printf("[ERROR] The number of seeds should be larger or equal to 1, aborting...\n");
		exit(EXIT_FAILURE);
	}
	chist->mtype = mtype;
	chist->k = k;
	chist->nseed = nseed;
	chist->xbin = xbin;
	chist->ybin = ybin;
	chist->ntot = xbin*ybin;
	chist->mbin = mbin;
	strncpy(chist->title, title, 79); //Maximum size of title.
	chist->dom = (Domain*)calloc(1, sizeof(Domain)); //Allocates the domain boundaries (xmin/xmax/ymin/ymax) without setting them.
	chist->colrule = (ColorRule*)calloc(1, sizeof(ColorRule));
	chist->data = (double*)calloc(chist->ntot, sizeof(double)); //The data should be initialized to zero everywhere.
	chist->mdata = (double*)calloc(chist->mbin, sizeof(double)); //The data should be initialized to zero everywhere.
}

/**
 * Deletes the content withing the given complex histogram, not the pointer "chist" itself.
 * Do not forget to invoke after use of "init_*" functions.
 */
void del_chistogram(Chistogram* chist) {
	free(chist->dom);
	free(chist->colrule);
	free(chist->data);
	free(chist->mdata);
}

/**
 * Returns the string corresponding to the given matrix type.
 */
char* matrixtype_to_string(MatrixType mtype) {//Private function
	switch (mtype) {
		case msmatrix:
			return "msmatrix";
		case normalized:
			return "normalized";
		default:
			return "none";
	}
}

/**
 * Converts the given string to the corresponding shape.
 * Returns 1 if the string is not recognized. 
 */
int string_to_matrixtype(const char* str, MatrixType* mtype) {//Private function
	if (strcmp(str, "msmatrix") == 0) {
		*mtype = msmatrix;
	}
	else if (strcmp(str, "normalized") == 0) {
		*mtype = normalized;
	}
	else {//Type not found.
		return 1;
	}
	return 0;
}

/**
 * Computes the first two statistical moments of the given list of complex numbers.
 */
void stat_moments(Chistogram* chist, int nz, dcomplex* z) {//Private function
	dcomplex avgz = 0.;
	int i;
	for (i = 0; i < nz; i++) {
		avgz += z[i];
	}
	avgz /= nz;
	chist->avgz = avgz;
	dcomplex dz;
	double varrez = 0., varimz = 0.;  //Variance in absolute value.
	for (i = 0; i < nz; i++) {
		dz = z[i] - avgz;
		varrez += creal(dz)*creal(dz);
		varimz += cimag(dz)*cimag(dz);
	}
	varrez /= nz;
	varimz /= nz;
	chist->sigmaz = sqrt(varrez) + I*sqrt(varimz); //Standard deviation.
	//printf("[INFO] Computed moments: avgz = %g%+gi, sigmaz = %g%+gi \n", creal(chist->avgz), cimag(chist->avgz), creal(chist->sigmaz), cimag(chist->sigmaz)); fflush(stdout);
}

/**
 * Initializes the bounds of the given complex histogram with a representative list of samples.
 * Does not take the values of the given list in account to feed the histogram.
 */
void setup_domain_chist(Chistogram* chist, int nz, dcomplex* z) {//Private function
	if (nz < 2) {
		fprintf(stderr, "[ERROR] The number of complex values initializing the histogram is too small, aborting...\n");
		exit(EXIT_FAILURE);
	}
	double x, xmin = creal(z[0]), xmax = creal(z[0]);
	double y, ymin = cimag(z[0]), ymax = cimag(z[0]);
	int i;
	for (i = 1; i < nz; i++) {//Find the extreme values of the given list
		x = creal(z[i]);
		y = cimag(z[i]);
		if (x < xmin)
			xmin = x;
		if (x > xmax)
			xmax = x;
		if (y < ymin)
			ymin = y;
		if (y > ymax)
			ymax = y;
	}
	if (!chist->dom_ok) {//If not set up yet, then initializes the domain bounds.
		chist->dom->xmin = xmin; //Extreme values may affect this. In this case, the domain should be setup explicitly by the user.
		chist->dom->xmax = xmax;
		chist->dom->ymin = ymin;
		chist->dom->ymax = ymax;
		chist->dom_ok = 1;  //Declare the domain as successfully set up.
		//printf("[INFO] The domain is set to xrange=%g:%g, yrange=%g:%g.\n", chist->xmin, chist->xmax, chist->ymin, chist->ymax);
	}
	if (!chist->mrange_ok) {//If not set up yet, then initializes the marginal distribution.
		chist->mmin = ymin;
		chist->mmax = ymax;
		chist->mrange_ok = 1;  //Declare the marginal bounds as successfully set up.
		//printf("[INFO] The marginal is set to xrange=%g:%g.\n", chist->mmin, chist->mmax);
	}
}

/**
 * Appends the given complex numbers "zlist" to the complex histogram "chist", after having set up the histogram domain.
 * Returns the number of complex numbers that lied outside of the histogram boundary, and which will not contribute.
 * @warning This function is the core of the histogram method. This algorithm MUST do the correct account according to the pixmap specification !
 */
int append_to_domain(Chistogram* chist, int nz, dcomplex* z) {//Private function
	if (!chist->dom_ok) {
		fprintf(stderr, "[ERROR] The domain boundaries are not set up, cannot append samples, aborting...\n");
		exit(EXIT_FAILURE);
	}
	int i, ix, iy, lost = 0;  //Number of missed numbers which will not contribute.
	for (i = 0; i < nz; i++) {//Loop on the values to append to the histogram.
		ix = (int)round((chist->xbin - 1)*(creal(z[i]) - chist->dom->xmin)/(chist->dom->xmax - chist->dom->xmin));
		iy = (int)round((chist->ybin - 1)*(chist->dom->ymax - cimag(z[i]))/(chist->dom->ymax - chist->dom->ymin));
		if (ix >= 0 && ix < chist->xbin && iy >= 0 && iy < chist->ybin) //If the value is in the domain, then add it.
			chist->data[ix + iy*chist->xbin] += 1.;
		else
			lost++;
	}
	return lost;
}

/**
 * Noramlizes the given domain histogram.
 */
void normalize_domain_chist(Chistogram* chist) {//Private function
	if (!chist->dom_ok) {
		fprintf(stderr, "[ERROR] The domain boundaries are not set up, cannot normalize, aborting...\n");
		exit(EXIT_FAILURE);
	}
	double total = 0.; //Total number of samples.
	int l;
	for (l = 0; l < chist->ntot; l++) {//Loop on the bins.
		total += chist->data[l];
	}
	double dx = (chist->dom->xmax - chist->dom->xmin)/(chist->xbin-1);
	double dy = (chist->dom->ymax - chist->dom->ymin)/(chist->ybin-1);
	for (l = 0; l < chist->ntot; l++) {//Normalize each bin.
		chist->data[l] /= dx*dy*total; //Each bin must be divided by its width.
	}
}

/**
 * Appends a list of samples to the marginal distribution of imaginary parts. Returns the number of lost samples.
 */
int append_to_marginal(Chistogram* chist, int nz, dcomplex* z) {//Private function
	if (!chist->mrange_ok) {
		fprintf(stderr, "[ERROR] The marginal boundaries are not set up, cannot append samples, aborting...\n");
		exit(EXIT_FAILURE);
	}
	double logx;
	double logxmin = log(chist->mmin);
	double logxmax = log(chist->mmax);
	int i, ix, lost = 0;
	for (i = 0; i < nz; i++) {//Loop on the samples to account for.
		logx = log(cimag(z[i])); //Already takes the logarithm here. This assumes positive imaginary part only.
		ix = (int)round((chist->mbin-1)*(logx - logxmin)/(logxmax - logxmin));
		if (0 <= ix && ix < chist->mbin) //Accounts for the sample only in the interval.
			chist->mdata[ix] += 1.;
		else //If not in the interval, then it is lost.
			lost++;
	}
	return lost;
}

/**
 * Normalize the given histogram to make the integral of the marginal being equal to 1. The logarithmic bin size is also accounted for.
 * This function must be called once immediately after finishing to append data to the histogram, and before any other operation.
 */
void normalize_marginal_chist(Chistogram* chist) {//Private function
	if (!chist->mrange_ok) {
		fprintf(stderr, "[ERROR] The marginal boundaries are not set up, cannot normalize, aborting...\n");
		exit(EXIT_FAILURE);
	}
	double total = 0.; //Total number of samples.
	int i;
	for (i = 0; i < chist->mbin; i++) {//Loop on the bins.
		total += chist->mdata[i];
	}
	double x, fac = log(chist->mmax/chist->mmin)/(chist->mbin-1);
	for (i = 0; i < chist->mbin; i++) {//Normalize each bin.
		x = chist->mmin*exp(fac*i); //Coordinate of the center of the i-th bin.
		chist->mdata[i] /= fac*x*total; //Each bin must be divided by its width.
	}
}

/**
 * Appends the list of eigenvalues to both the domain histogram and the marginal histogram.
 */
int fill_chist(Chistogram* chist, int nz, dcomplex* z) {
	int lost = 0; //Total number of lost eigenvalues.
	stat_moments(chist, nz, z); //Computes the statistical moments.
	setup_domain_chist(chist, nz, z);  //Setup the histogram domain.
	lost += append_to_domain(chist, nz, z);   //Appends to the rectangular domain.
	lost += append_to_marginal(chist, nz, z); //Appends to the marginal distribution of imaginary parts.
	//normalize_domain_chist(chist);    //Normalizes the probability density of the rectangular histogram (optional).
	normalize_marginal_chist(chist);  //Normalizes the probability density in logarithmic coordinates (required).
	return lost;
}

/**
 * Initializes the minimum and maximum bin values in the color rule for later use.
 */
void find_min_max_chist(Chistogram* chist, int verbose) {//Private function
	int l, hmin = chist->data[0], hmax = chist->data[0];
	for (l = 1; l < chist->ntot; l++) {
		if (hmin > chist->data[l])
			hmin = chist->data[l];
		if (hmax < chist->data[l])
			hmax = chist->data[l];
	}
	if (hmax == hmin) {//Protects against zero division.
		if (verbose)
			printf("[WARN] The histogram does not contain useful data (hmin=hmax=%d)...\n", hmin);
		hmax += 1.;
	}
	chist->colrule->hmin = hmin;
	chist->colrule->hmax = hmax;
	//printf("[INFO] Setting color rule with hmin=%d, hmax=%d\n", hmin, hmax);
}

/**
 * Prints a short summary of the parameters used in the complex histogram for the standard output.
 */
void print_param_chist(Chistogram* chist) {
	printf("[INFO] Complex histogram entitled '%s' for matrix type '%s' with k=%g%+gi, size %dx%d and %d seeds.\n",
		chist->title, matrixtype_to_string(chist->mtype), creal(chist->k), cimag(chist->k), chist->xbin, chist->ybin, chist->nseed);
}

/**
 * Exports the complex hisogram to image PPM file "fname.ppm", then to PNG file "fname.png".
 * Note that the filename "fname" should not have any extension.
 */
void export_image_chistogram(Chistogram* chist, int id, const char* fname) {//Private function
	size_t imgfnamelen = strlen(fname) + (size_t)floor(1+log10(1+abs(id))) + 2;
	char imgfname[imgfnamelen];
	sprintf(imgfname, "%s_%d", fname, id);
	find_min_max_chist(chist, 1); //Find min/max values (1=verbose mode).
	char ppmfile[strlen(imgfname)+6];
	sprintf(ppmfile, "%s.ppm", imgfname);
	FILE* outfp = fopen(ppmfile, "wb");
	fprintf(outfp, "P6\n%d %d\n255\n", chist->xbin, chist->ybin);  //Header with 'P6' for binary PPM file format.
	Pixel pix;
	int l;
	for (l = 0; l < chist->ntot; l++) {//Assume logarithmic scale.
		set_color(chist->colrule, chist->data[l], &pix);
		fwrite(&pix, 1, sizeof(Pixel), outfp);
	}
	fclose(outfp);
	convert_ppm_to_png(imgfname);
}

/**
 * Creates the PNG image "fname_id.png" and appends a "tikzpicture" environment importing "fname_id.png" to the given TikZ file "fname.tikz".
 * The filename "fname" should not have any extension.
 */
void export_tikz_chistogram(Chistogram* chist, Medium* med, int id, const char* fname) {
	if (!chist->dom_ok || !chist->mrange_ok) {
		fprintf(stderr, "[ERROR] The boundaries are not set up, cannot plot, aborting...\n");
		exit(EXIT_FAILURE);
	}
	export_image_chistogram(chist, id, fname);
	char tikzfname[strlen(fname)+6], varname[20];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a");
	double dx = (chist->dom->xmax - chist->dom->xmin)/(chist->xbin - 1);
	double dy = (chist->dom->ymax - chist->dom->ymin)/(chist->ybin - 1);
	double txmin = chist->dom->xmin - dx/2; //True TikZ min/max ranges.
	double txmax = chist->dom->xmax + dx/2; //Indeed, the bounds of the complex map consider the center of each pixels.
	double tymin = chist->dom->ymin - dy/2;
	double tymax = chist->dom->ymax + dy/2;
	fprintf(tikzfp, "\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n\tname={plane},\n");
	fprintf(tikzfp, "\txmin=%g, xmax=%g,\n\tymin=%g, ymax=%g,\n", txmin, txmax, tymin, tymax);
	if (chist->mtype == msmatrix) {
		sprintf(varname, "\\mu");
		fprintf(tikzfp, "\ttitle={%dD %s with $N=%d$ and %s\\\\ for $k=(%g%+g\\,{\\rm i})~\\varsigma^{-1}$ with %d seeds\\\\ map $%d\\times %d$ (%.4g s)},\n\ttitle style={align=center},\n",
			med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, creal(chist->k), cimag(chist->k), chist->nseed, chist->xbin, chist->ybin, chist->realtime);
	}
	else {//Normalized MS matrix.
		sprintf(varname, "\\nu");
		fprintf(tikzfp, "\ttitle={%dD %s with $N=%d$\\\\ for $k=(%g%+g\\,{\\rm i})~\\varsigma^{-1}$ with %d seeds\\\\ map $%d\\times %d$ (%.4g s)},\n\ttitle style={align=center},\n",
			med->d, shape_to_string(med->shape), med->n, creal(chist->k), cimag(chist->k), chist->nseed, chist->xbin, chist->ybin, chist->realtime);
	}
	fprintf(tikzfp, "\txlabel={$\\Re%s$},\n\tylabel={$\\Im%s$},\n", varname, varname);
	fprintf(tikzfp, "\taxis equal image,\n\tyticklabel style={rotate=90},\n\taxis on top,\n\tlegend pos={north west},\n");
	char* colstr = tikz_colorbar_v2(chist->colrule);
	fprintf(tikzfp, "%s", colstr);
	free(colstr);
	fprintf(tikzfp, "]%%\n\\addplot[forget plot] graphics[xmin=%g, xmax=%g, ymin=%g, ymax=%g] {%s_%d.png};\n", txmin, txmax, tymin, tymax, fname, id);
	//fprintf(tikzfp, "\\addplot[black, domain=0:360, samples=32, smooth] ({%g+%g*cos(x)},{%g+%g*sin(x)});\n", creal(chist->avgz), chist->sigmaz, cimag(chist->avgz), chist->sigmaz);
	fprintf(tikzfp, "\\fill[black] (axis cs:%g,%g) circle[radius=1.2pt]; %%Average point\n", creal(chist->avgz), cimag(chist->avgz));
	int i, nsamp = 150;
	double r, rmin = 0.001, rmax = 1.; //Radial coordinate (about rmax=1 is generally enough).
	double dr = (rmax-rmin)/(nsamp-1);
	dcomplex mu;
	fprintf(tikzfp, "\\addplot[black, densely dotted, thick] coordinates {%%Green spiral %s+ (symmetric two-atom states)\n\t", varname);
	for (i = 0; i < nsamp; i++) {//Draw the Green spiral "mu+"
		r = rmin + i*dr;
		if (chist->mtype == msmatrix)
			mu = invf(med, chist->k) - free_green(med->d, chist->k, r);
		else //Normalized MS matrix
			mu = I - free_green(med->d, chist->k, r)/free_green_imag(med->d, chist->k, 0.);
		fprintf(tikzfp, "(%.16g, %.16g) ", creal(mu), cimag(mu));
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{$%s_+$}\n", varname);
	fprintf(tikzfp, "\\addplot[black, densely dashed, thick] coordinates {%%Green spiral %s- (antisymmetric two-atom states)\n\t", varname);
	for (i = 0; i < nsamp; i++) {//Draw the Green spiral "mu-"
		r = rmin + i*dr;
		if (chist->mtype == msmatrix)
			mu = invf(med, chist->k) + free_green(med->d, chist->k, r);
		else //Normalized MS matrix
			mu = I + free_green(med->d, chist->k, r)/free_green_imag(med->d, chist->k, 0.);
		fprintf(tikzfp, "(%.16g, %.16g) ", creal(mu), cimag(mu));
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{$%s_-$}\n", varname);
	fprintf(tikzfp, "\\end{axis}%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\tname={marginal},\n\tat={(plane.below south west)},\n\tanchor={above north west},\n\tyshift=-5pt,\n");
	fprintf(tikzfp, "\txmin=%g, xmax=%g,\n", chist->mmin, chist->mmax); //\txmode=log, ymode=log,\n
	fprintf(tikzfp, "\ttitle={(b) Marginal distribution of $\\Im%s$},\n\ttitle style={align=center},\n", varname);
	fprintf(tikzfp, "\txlabel={$\\Im%s$},\n\tylabel={$p(\\Im%s)$},\n", varname, varname);
	fprintf(tikzfp, "\tyticklabel style={rotate=90},\n\taxis on top=true,\n\tlegend pos={outer north east},\n]%%\n");
	fprintf(tikzfp, "\\addplot[blue, thick] coordinates {%%Marginal histogram Im(%s), avgx=%g, sigmax=%g\n\t", varname, cimag(chist->avgz), cimag(chist->sigmaz));
	double x, fac = log(chist->mmax/chist->mmin)/(chist->mbin-1);
	for (i = 0; i < chist->mbin; i++) {//Loop on the samples to plot.
		x = chist->mmin*exp(fac*i); //Coordinate of the center of the i-th bin.
		fprintf(tikzfp, "(%g, %g) ", x, chist->mdata[i]);
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{Histogram}\n");
	double sigmax = cimag(chist->sigmaz); //Exact, but more accurate to compute this from Tr(G⁺G).
	//double sigmax_expc = 0.5*sqrt(med->n-1)*(3./(8*PI*radius_of_ball(med)))/free_green_imag(med->d, chist->k, 0.); //Std deviation of "nu" in 3D.
	//printf("[INFO] Theoretical sigma=%g, exact sigma=%g \n", sigmax_expc, sigmax);
	double avgx = cimag(chist->avgz);
	double xp = pow(avgx + sigmax, 2)/avgx; //Estimates the Marchenko-Pastur parameters from the theoretical mean and variance.
	double xm = pow(avgx - sigmax, 2)/avgx;
	double px, cn = 0.5*PI*pow(sqrt(xp)-sqrt(xm), 2); //Normalization coefficient located in the denominator.
	fprintf(tikzfp, "\\addplot[black, densely dotted, thick] coordinates {%%Marchenko-Pastur distribution (xm=%g, xp=%g)\n\t", xm, xp);
	fac = PI/(nsamp-1);
	for (i = 1; i <= nsamp-2; i++) {//Loop on the samples to plot.
		x = ((xp + xm) - (xp - xm)*cos(fac*i))/2;  //Position of the i-th interval using more appropriate Chebyshev samples.
		px = sqrt((xp-x)*(x-xm))/(cn*x); //Marchenko-Pastur density itself.
		fprintf(tikzfp, "(%g, %g) ", x, px);
	}
	fprintf(tikzfp, "\n}; \\addlegendentry{Marchenko-Pastur}\n");
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	fclose(tikzfp);
}

/**
 * Saves the complex histogram data, in addition to other metadata, to the file "fname".
 */
void save_chistogram(Chistogram* chist, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[muplane]\ntitle=%s\n", chist->title);
	fprintf(fp, "mtype=%s\n", matrixtype_to_string(chist->mtype));
	fprintf(fp, "k=%.16g%+.16gi\n", creal(chist->k), cimag(chist->k));
	fprintf(fp, "nseed=%d\n", chist->nseed);
	fprintf(fp, "xbin=%d\n", chist->xbin);
	fprintf(fp, "ybin=%d\n", chist->ybin);
	save_domain(chist->dom, fp);
	fprintf(fp, "avgz=%.16g%+.16gi\n", creal(chist->avgz), cimag(chist->avgz));
	fprintf(fp, "sigmaz=%.16g%+.16gi\n", creal(chist->sigmaz), cimag(chist->sigmaz));
	fprintf(fp, "realtime=%g\n", chist->realtime);
	save_colorrule(chist->colrule, "color", fp);
	fprintf(fp, "data=");
	int l;
	for (l = 0; l < chist->ntot; l++) {
		fprintf(fp, "%.16g ", chist->data[l]);
	}
	fprintf(fp, "\nmbin=%d\n", chist->mbin);
	fprintf(fp, "mrange=%.16g:%.16g\n", chist->mmin, chist->mmax);
	fprintf(fp, "mdata=");
	for (l = 0; l < chist->mbin; l++) {
		fprintf(fp, "%.16g ", chist->mdata[l]);
	}
	fprintf(fp, "\n\n");
	fclose(fp);
}

/**
 * Parses the complex histogram "chist" with the arguments of "args".
 * The "chist" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_chistogram(Chistogram* chist, int narg, char** args) {
	int nargmin = 6; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	char *value, *title, c;
	MatrixType mtype;
	double re, im, mmin, mmax;
	int xbin, ybin, mbin, nseed, ncount;
	dcomplex k;
	value = get_value(narg, args, "k");
	if (sscanf(value, "%lg%lg%c", &re, &im, &c) == 3 && c == 'i') {
		k = re + im*I;
	}
	else {
		printf("[ERROR] Invalid complex wave number '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "nseed");
	if (sscanf(value, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "mtype");
	if (string_to_matrixtype(value, &mtype)) {
		printf("[ERROR] Invalid matrix type '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "xbin");
	if (sscanf(value, "%d", &xbin) != 1) {
		printf("[ERROR] Invalid number of horizontal bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "ybin");
	if (sscanf(value, "%d", &ybin) != 1) {
		printf("[ERROR] Invalid number of vertical bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "mbin");
	if (sscanf(value, "%d", &mbin) != 1) {
		printf("[ERROR] Invalid number of marginal bins '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	title = get_value(narg, args, "title");
	init_chistogram(chist, mtype, k, nseed, xbin, ybin, mbin, title); //Call initialization and allocation function.
	chist->dom_ok = !parse_domain(chist->dom, narg, args); //parse_domain() returns 1 on failure and 0 on success.
	if (!chist->dom_ok) {
		printf("[WARN] The domain boundaries are not set up yet...\n");
	}
	value = get_value(narg, args, "mrange");
	if (value[0]) {//If given field "mrange"
		if (sscanf(value, "%lg:%lg", &mmin, &mmax) != 2) {
			printf("[ERROR] Invalid marginal range '%s', aborting...\n", value);
			exit(EXIT_FAILURE);
		}
		chist->mmin = mmin;
		chist->mmax = mmax;
		chist->mrange_ok = 1;
	}
	else {//No "range" found.
		chist->mrange_ok = 0;
	}
	if (sscanf(get_value(narg, args, "avgz"), "%lg%lg%c", &re, &im, &c) == 3 && c == 'i') {
		chist->avgz = re + im*I;
	}
	if (sscanf(get_value(narg, args, "sigmaz"), "%lg%lg%c", &re, &im, &c) == 3 && c == 'i') {
		chist->sigmaz = re + im*I;
	}
	sscanf(get_value(narg, args, "realtime"), "%lg", &(chist->realtime));
	value = get_value(narg, args, "color"); //Now parses the colorscheme.
	parse_colorrule(chist->colrule, value);
	value = get_value(narg, args, "data");  //Now parses the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != chist->ntot) {
			fprintf(stderr, "[ERROR] Histogram data corruption. Expected %d data points but found %d, aborting...\n", chist->ntot, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(chist->ntot, chist->data, value);
	}
	value = get_value(narg, args, "mdata");  //Now parses the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != chist->mbin) {
			fprintf(stderr, "[ERROR] Histogram data corruption. Expected %d data points but found %d, aborting...\n", chist->mbin, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(chist->mbin, chist->mdata, value);
	}
}
