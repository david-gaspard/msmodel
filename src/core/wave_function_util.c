/****
 * @date Created on 2021-11-27 at 18:31:12 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to manipulate and plot the wave function structure.
 ***/
#include <stdlib.h>               /* Standard Library for Memory Allocation */
#include <stdio.h>                /* Standard Library for Input and Output */
#include <string.h>               /* Standard Library for String Manipulation */
#include <math.h>                 /* Standard Library for Mathematical Functions */
#include "common_util.h"          /* Import the General Purpose Functions */
#include "wave_function_type.h"   /* Import the Wave Function Type */
#include "domain_util.h"          /* Import the 2D Domain Utilities */
#include "incident_wave_util.h"   /* Import the Incident Wave Utilities */
#include "real_vector_util.h"     /* Import the Complex Data Parsing functions */
#include "medium_util.h"          /* Import the Medium Management Utilities */
#include "color_rule_util.h"      /* Import the Color Rule Utilities */

/**
 * Initialization function of an empty complex map with the given parameters.
 * The content of the complex map must be deleted with the corresponding del_*() function.
 */
void init_wfun(WaveFunction* wfun, int nx, int ny, Domain* dom, int nseed, uint64_t iseed, IncidentWave* iwave) {
	if (nx*ny > 10000000) {//Beyond 10 MegaPixels is a bit too large...
		printf("[ERROR] The number of samples is too large (nx=%d, ny=%d), aborting...\n", nx, ny);
		exit(EXIT_FAILURE);
	}
	if (nx < 2 || ny < 2) {
		printf("[ERROR] The number of samples is too small (nx=%d, ny=%d), aborting...\n", nx, ny);
		exit(EXIT_FAILURE);
	}
	if (nseed <= 0) {
		printf("[ERROR] Invalid number of seeds (nseed=%d), aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	wfun->nx = nx;
	wfun->ny = ny;
	wfun->dx = (dom->xmax - dom->xmin)/(nx - 1);
	wfun->dy = (dom->ymax - dom->ymin)/(ny - 1);
	wfun->dom = (Domain*)calloc(1, sizeof(Domain));
	copy_domain(dom, wfun->dom);
	wfun->iwave = (IncidentWave*)calloc(1, sizeof(IncidentWave));
	copy_iwave(iwave, wfun->iwave);
	wfun->nseed = nseed;
	wfun->iseed = iseed;
	wfun->realtime = 0.; //Computation time is zero, since empty.
	int ntot = nx*ny;
	wfun->function = (double*)calloc(ntot, sizeof(double));
	wfun->density = (double*)calloc(ntot, sizeof(double));
	wfun->colrule1 = (ColorRule*)calloc(1, sizeof(ColorRule));
	wfun->colrule2 = (ColorRule*)calloc(1, sizeof(ColorRule));
}

/**
 * Deletes the content (allocated pointers) within the given wave function, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_wfun(WaveFunction* wfun) {
	free(wfun->dom);
	free(wfun->iwave);
	free(wfun->function);
	free(wfun->density);
	free(wfun->colrule1);
	free(wfun->colrule2);
}

/**
 * Returns the 2D position corresponding to the given index "l" in the range {0, ..., nx*ny-1}.
 * Note that the other directions in "pos" are unchanged, and typically are left at zero.
 * The samples are assumed to be stored in row-major format in wfun->function and wfun->density.
 */
void get_pos_wfun(WaveFunction* wfun, int l, double* pos) {
	int ix = l%wfun->nx;
	int iy = l/wfun->nx;
	pos[0] = (wfun->dom->xmin + ix*wfun->dx);
	pos[1] = (wfun->dom->ymax - iy*wfun->dy);
}

/**
 * Divides the arrays of the wave function by "div".
 */
void divide_wfun(WaveFunction* wfun, double div) {
	int l, ntot = wfun->nx*wfun->ny;
	for (l = 0; l < ntot; l++) {//Loop on the pixels.
		wfun->function[l] /= div;
		wfun->density[l] /= div;
	}
}

/**
 * Prints a short summary of the parameters used in the wave function for the standard output.
 */
void print_param_wfun(WaveFunction* wfun) {
	printf("[INFO] Wave function of size %dx%d in region x=%g:%g y=%g:%g for %d seeds with source %s.\n",
		wfun->nx, wfun->ny, wfun->dom->xmin, wfun->dom->xmax, wfun->dom->ymin, wfun->dom->ymax, wfun->nseed, wfun->iwave->desc_ascii);
}

/**
 * Sets the minimum and maximum values of real part found in the wave function.
 * The "hmin" and "hmax" values are stored in the color rules for later use.
 */
void set_min_max_wfun(WaveFunction* wfun, int verbose) {
	double h1min = wfun->function[0], h1max = wfun->function[0]; //Initialize the min/max values for the function map.
	double h2min = wfun->density[0],  h2max = wfun->density[0];  //Initialize the min/max values for the density map.
	int l, ntot = wfun->nx*wfun->ny;
	for (l = 1; l < ntot; l++) {
		if (wfun->function[l] < h1min) {
			h1min = wfun->function[l];
		}
		if (wfun->function[l] > h1max) {
			h1max = wfun->function[l];
		}
		if (wfun->density[l] < h2min) {
			h2min = wfun->density[l];
		}
		if (wfun->density[l] > h2max) {
			h2max = wfun->density[l];
		}
	}
	if (h1min == h1max) {//Protects against zero division.
		if (verbose) {
			printf("[WARN] The wave function map does not contain useful data (h1min=h1max=%g)...\n", h1min);
		}
		h1max += 1.;
	}
	if (h2min == h2max) {//Protects against zero division.
		if (verbose) {
			printf("[WARN] The wave function density does not contain useful data (h2min=h2max=%g)...\n", h2min);
		}
		h2max += 1.;
	}
	wfun->colrule1->hmin = h1min; //Stores the minimum and maximum real parts in the respective color rule of the two maps.
	wfun->colrule1->hmax = h1max;
	wfun->colrule2->hmin = h2min;
	wfun->colrule2->hmax = h2max;
}

/**
 * Plots both the real part and the density arrays of the wave function "wfun" of the given wave function "wfun" to a PPM file (Portable PixMap).
 * This function then attempts to convert the PPM file to a PNG file if the appropriate command "pnmtopng" exists.
 * Note that the filename "fname" should not have any extension.
 */
void export_image_wfun(WaveFunction* wfun, int id, const char* fname) {
	set_min_max_wfun(wfun, 0); //Find maximum values (0=quiet mode).
	size_t imgfnamelen = strlen(fname) + (size_t)floor(1+log10(1+abs(id))) + 3;
	char imgfname1[imgfnamelen], imgfname2[imgfnamelen];
	sprintf(imgfname1, "%s_%d_1", fname, id);
	sprintf(imgfname2, "%s_%d_2", fname, id);
	char ppmfile1[strlen(imgfname1)+6], ppmfile2[strlen(imgfname2)+6];
	sprintf(ppmfile1, "%s.ppm", imgfname1);
	sprintf(ppmfile2, "%s.ppm", imgfname2);
	FILE* fp1 = fopen(ppmfile1, "wb");  //'wb' for binary writing.
	FILE* fp2 = fopen(ppmfile2, "wb");  //'wb' for binary writing.
	fprintf(fp1, "P6\n%d %d\n255\n", wfun->nx, wfun->ny);  //Header with 'P6' for binary PPM file format.
	fprintf(fp2, "P6\n%d %d\n255\n", wfun->nx, wfun->ny);  //Header with 'P6' for binary PPM file format.
	Pixel pix1, pix2;
	int l, ntot = wfun->nx*wfun->ny;
	for (l = 0; l < ntot; l++) {//Loop on pixels to write.
		set_color(wfun->colrule1, wfun->function[l], &pix1); //This assumes "hmin" and "hmax" are already initialized (see set_min_max_wfun).
		set_color(wfun->colrule2, wfun->density[l],  &pix2); //This assumes "hmin" and "hmax" are already initialized (see set_min_max_wfun).
		fwrite(&pix1, 1, sizeof(Pixel), fp1);
		fwrite(&pix2, 1, sizeof(Pixel), fp2);
	}
	fclose(fp1);
	fclose(fp2);
	convert_ppm_to_png(imgfname1);
	convert_ppm_to_png(imgfname2);
}

/**
 * Appends a TikZ picture environment importing the image "fname.png" to a file "fname".
 * Add the metadata in "wfun" to the plot, and a color legend.
 * Note that the filename "fname" should not have any extension.
 */
void export_tikz_wfun(WaveFunction* wfun, Medium* med, int id, const char* fname) {
	export_image_wfun(wfun, id, fname);
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	double txmin = wfun->dom->xmin - wfun->dx/2; //True TikZ min/max ranges.
	double txmax = wfun->dom->xmax + wfun->dx/2; //Indeed, the bounds of the complex map consider the center of each pixels.
	double tymin = wfun->dom->ymin - wfun->dy/2;
	double tymax = wfun->dom->ymax + wfun->dy/2;
	fprintf(tikzfp, "\\par%%\n\\begin{tikzpicture}[%%\n\t/pgfplots/width=360pt,\n\t/pgfplots/xmin=%g,\n\t/pgfplots/xmax=%g,\n\t/pgfplots/ymin=%g,\n\t/pgfplots/ymax=%g,\n\t/pgfplots/xlabel={$x~(\\varsigma)$},\n\t/pgfplots/ylabel={$y~(\\varsigma)$},\n\t/pgfplots/title style={align=center},\n\t/pgfplots/enlargelimits=false,\n\t/pgfplots/yticklabel style={rotate=90},\n\t/pgfplots/axis on top=true,\n]%%\n",
		txmin, txmax, tymin, tymax);
	fprintf(tikzfp, "\\begin{axis}[%%\n\tname={wfun-function},\n\ttitle={%dD %s for $N=%d$ and %s\\\\ $\\Re\\langle\\psi(\\mathbf{r})\\rangle$, source %s\\\\ with %d seeds (map $%d\\times %d$, %.1f s)},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, wfun->iwave->desc_latex, wfun->nseed, wfun->nx, wfun->ny, wfun->realtime);
	fprintf(tikzfp, "\taxis equal image,\n");
	fprintf(tikzfp, "%s", tikz_colorbar_v1(wfun->colrule1)); //Colorbar v1=equally spaced ticks, v2=ticks are powers of ten.
	fprintf(tikzfp, "]%%\n\\addplot[forget plot] graphics[xmin=%g, xmax=%g, ymin=%g, ymax=%g] {%s_%d_1.png};\n", txmin, txmax, tymin, tymax, fname, id);
	fprintf(tikzfp, "\\addplot[black, opacity=0.4, domain=0:360, samples=64, smooth] ({%g*cos(x)}, {%g*sin(x)});\n", radius_of_ball(med), radius_of_ball(med));
	fprintf(tikzfp, "\\end{axis}%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\tname={wfun-density},\n\tat={(wfun-function.below south west)},\n\tanchor={above north west},\n\tyshift=-10pt,\n\ttitle={%dD %s for $N=%d$ and %s\\\\ $\\langle|\\psi(\\mathbf{r})|^2\\rangle$, source %s\\\\ with %d seeds (map $%d\\times %d$, %.1f s)},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, wfun->iwave->desc_latex, wfun->nseed, wfun->nx, wfun->ny, wfun->realtime);
	fprintf(tikzfp, "\taxis equal image,\n");
	fprintf(tikzfp, "%s", tikz_colorbar_v1(wfun->colrule2)); //Colorbar v1=equally spaced ticks, v2=ticks are powers of ten.
	fprintf(tikzfp, "]%%\n\\addplot[forget plot] graphics[xmin=%g, xmax=%g, ymin=%g, ymax=%g] {%s_%d_2.png};\n", txmin, txmax, tymin, tymax, fname, id);
	if (med->shape == ball) {//Show the shape of the medium.
		double radius = radius_of_ball(med); //Pre-computes the radius of the ball.
		fprintf(tikzfp, "\\addplot[black, opacity=0.4, domain=0:360, samples=64, smooth] ({%g*cos(x)}, {%g*sin(x)});\n", radius, radius);
	}
	if (med->d == 2 && wfun->nseed == 1) {//Show the points only in 2D for a single geometry.
		fill_medium(med, wfun->iseed); //Ensures the medium is filled.
		fprintf(tikzfp, "\\addplot[black, mark=*, only marks, mark options={black, scale=0.15, fill opacity=0.7, draw opacity=0}] coordinates {%% Scatterer positions.");
		int i, np = 0;  //Counts the number of plotted points.
		for (i = 0; i < med->n; i++) {//Loop on atoms.
			if (is_point_in_domain(wfun->dom, med->pos+2*i)) {//If the point is in the domain, then show it.
				if (np%20 == 0) {//Reduce line length.
					fprintf(tikzfp, "\n\t");
				}
				fprintf(tikzfp, "(%.12g,%.12g) ", med->pos[2*i], med->pos[2*i+1]);
				np++;
			}
		}
		fprintf(tikzfp, "\n}; %% Plotted %d/%d points.\n", np, med->n);
	}
	fprintf(tikzfp, "\\end{axis}%%\n\\end{tikzpicture}%%\n");
	fclose(tikzfp);
}

/**
 * Saves the wave function data, in addition to other metadata, to the file titled "fname".
 */
void save_wfun(WaveFunction* wfun, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[wavefunction]\n");
	fprintf(fp, "nseed=%d\n", wfun->nseed);
	fprintf(fp, "iseed=%lu\n", wfun->iseed);
	fprintf(fp, "nx=%d\n", wfun->nx);
	fprintf(fp, "ny=%d\n", wfun->ny);
	save_domain(wfun->dom, fp);
	save_iwave(wfun->iwave, fp);
	fprintf(fp, "realtime=%g\n", wfun->realtime); //Total physical duration of the computation of the complex map (in seconds).
	save_colorrule(wfun->colrule1, "color1", fp);
	save_colorrule(wfun->colrule2, "color2", fp);
	int ntot = wfun->nx*wfun->ny;
	save_real_data(ntot, wfun->function, "function", fp);
	save_real_data(ntot, wfun->density, "density", fp);
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Parses the wave function "wfun" with the arguments of "args".
 * The "wfun" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_wfun(WaveFunction* wfun, int narg, char** args) {
	int nargmin = 8; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	int nx = 0;  //Default value.
	int ny = 0;  //Default value.
	int nseed;
	uint64_t iseed;
	char *value;
	value = get_value(narg, args, "nx");
	if (sscanf(value, "%d", &nx) != 1) {
		printf("[ERROR] Invalid number of sample nx='%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "ny");
	if (sscanf(value, "%d", &ny) != 1) {
		printf("[ERROR] Invalid number of sample ny='%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	Domain dom; //Parse the complex domain.
	if (parse_domain(&dom, narg, args)) {//If the domain parsing failed
		printf("[ERROR] Invalid domain, aborting...\n");
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
	init_wfun(wfun, nx, ny, &dom, nseed, iseed, &iwave);
	sscanf(get_value(narg, args, "realtime"), "%lg", &(wfun->realtime)); //Set realtime after initialization.
	parse_colorrule(wfun->colrule1, get_value(narg, args, "color1")); //Now parsing the first color scheme (function).
	parse_colorrule(wfun->colrule2, get_value(narg, args, "color2")); //Now parsing the second color scheme (density).
	int ncount, ncount_expc = nx*ny; //Expected number of real points.
	value = get_value(narg, args, "function");  //Now parsing the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != ncount_expc) {
			fprintf(stderr, "[ERROR] Wave function data corruption. Expected %d data points but found %d, aborting...\n", ncount_expc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ncount, wfun->function, value);
	}
	value = get_value(narg, args, "density");  //Now parsing the list of data.
	if (value[0]) {//If existing data field, then parses it.
		ncount = count_real_data(value);
		if (ncount != ncount_expc) {
			fprintf(stderr, "[ERROR] Wave function data corruption. Expected %d data points but found %d, aborting...\n", ncount_expc, ncount);
			exit(EXIT_FAILURE);
		}
		parse_real_data(ncount, wfun->density, value);
	}
}
