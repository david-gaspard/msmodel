/****
 * @date Created on 2020-08-11 at 18:09 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the functions to create and manipulate complex maps.
 ***/
#include <stdlib.h>                  /* Standard Library for Memory Allocation */
#include <stdio.h>                   /* Standard Library for Input and Output */
#include <string.h>                  /* Standard Library for String Manipulation */
#include <math.h>                    /* Standard Library for Mathematical Functions */
#include "common_util.h"             /* Import the General Purpose Functions */
#include "complex_map_type.h"        /* Import the Complex Map Type */
#include "complex_data_util.h"       /* Import the Complex Data Parsing functions */
#include "medium_util.h"             /* Import the Medium Management Utilities */
#include "color_rule_util.h"         /* Import the Color Rule Utilities */
#include "square_well_resonances.h"  /* Import the Routine to Find the Resonances in the Square Well Approximation */

/**
 * Returns the position in the z-plane corresponding to the given index "l" in the range {0, ..., nx*ny-1}.
 * The samples are assumed to be stored in row-major format.
 */
dcomplex get_arg_cmap(ComplexMap* cmap, int l) {
	int ix = l%cmap->nx;
	int iy = l/cmap->nx;
	return (cmap->xmin + ix*cmap->dx) + I*(cmap->ymax - iy*cmap->dy);
}

/**
 * Returns the element of the complex map at the given index (ix, iy). The origin is in the top left corner.
 * The index ix runs horizontally from left to right, and iy runs vertically from up to bottom.
 * Do not check for bounds.
 */
dcomplex get_value_cmap(ComplexMap* cmap, int ix, int iy) {
	return cmap->data[iy*cmap->nx + ix];
}

/**
 * Assigns the given value "val" at the given position (ix, iy) in the complex map. Do not make bound checking.
 */
void set_value_cmap(ComplexMap* cmap, int ix, int iy, dcomplex val) {
	cmap->data[iy*cmap->nx + ix] = val;
}

/**
 * Initialization function of an empty complex map with the given parameters.
 * The content of the complex map must be deleted with the corresponding del_*() function.
 */
void init_complexmap(ComplexMap* cmap, int nx, int ny, double xmin, double xmax, double ymin, double ymax, int nseed, const char* title) {
	if (xmax <= xmin) {
		printf("[ERROR] xmin=%f should be smaller than xmax=%f, aborting...\n", xmin, xmax);
		exit(EXIT_FAILURE);
	}
	if (ymax <= ymin) {
		printf("[ERROR] ymin=%f should be smaller than ymax=%f, aborting...\n", ymin, ymax);
		exit(EXIT_FAILURE);
	}
	if (nx*ny > 10000000) {//Beyond 10 MegaPixels is a bit too large...
		printf("[ERROR] The number of samples is too large (nx=%d, ny=%d), aborting...\n", nx, ny);
		exit(EXIT_FAILURE);
	}
	if ((nx < 2) != (ny < 2)) {//If one of the number of samples is invalid, then try to infer it from the other one.
		if (nx < 2) {//Infer nx from ny.
			nx = (int)round(1.0 + (ny-1)*(xmax-xmin)/(ymax-ymin));
		}
		else {//Infer ny from nx.
			ny = (int)round(1.0 + (nx-1)*(ymax-ymin)/(xmax-xmin));
		}
		printf("[INFO] Inferred map size %dx%d totalling %d samples...\n", nx, ny, nx*ny);
	}
	if (nx < 2 || ny < 2) {
		printf("[ERROR] The number of samples is too small (nx=%d, ny=%d), aborting...\n", nx, ny);
		exit(EXIT_FAILURE);
	}
	if (nseed <= 0) {
		printf("[ERROR] Invalid number of seeds (nseed=%d), aborting...\n", nseed);
		exit(EXIT_FAILURE);
	}
	cmap->nx = nx;
	cmap->ny = ny;
	cmap->ntot = nx*ny;
	cmap->xmin = xmin;
	cmap->xmax = xmax;
	cmap->ymin = ymin;
	cmap->ymax = ymax;
	cmap->dx = (xmax - xmin)/(nx - 1);
	cmap->dy = (ymax - ymin)/(ny - 1);
	cmap->nseed = nseed;
	strncpy(cmap->title, title, 79); //Maximum size of the title.
	cmap->realtime = 0.0; //Computation time is zero, since empty.
	cmap->data = (dcomplex*)calloc(cmap->ntot, sizeof(dcomplex));
	cmap->get_arg = (dcomplex(*)(void*, int))get_arg_cmap;
	cmap->colrule = (ColorRule*)calloc(1, sizeof(ColorRule));
}

/**
 * Deletes the content (allocated pointers) within the given complex map, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_complexmap(ComplexMap* cmap) {
	free(cmap->colrule);
	free(cmap->data);
}

/**
 * Sets the minimum and maximum values of real part found in the complex map.
 * The "hmin" and "hmax" values are stored in the color rule for later use.
 */
void find_min_max_cmap(ComplexMap* cmap, int verbose) {
	double hmin = creal(cmap->data[0]); //Initialize the minimum real part.
	double hmax = creal(cmap->data[0]); //Initialize the maximum real part.
	int l;
	for (l = 1; l < cmap->ntot; l++) {
		if (creal(cmap->data[l]) < hmin) {
			hmin = creal(cmap->data[l]);
		}
		if (creal(cmap->data[l]) > hmax) {
			hmax = creal(cmap->data[l]);
		}
	}
	if (hmin == hmax) {//Protects against zero division.
		if (verbose) {
			printf("[WARN] The complex map does not contain useful data (hmin=hmax=%g)...\n", hmin);
		}
		hmax += 1.0;
	}
	cmap->colrule->hmin = hmin; //Stores the minimum and maximum real parts in the color rule.
	cmap->colrule->hmax = hmax;
}

/**
 * Computes the value of the map at the given point "z" using bilinear interpolation.
 */
dcomplex interpolate_bilinear(ComplexMap* cmap, dcomplex z) {
	double wx = (creal(z)-cmap->xmin)/cmap->dx;  //Continuous indices.
	double wy = (cmap->ymax-cimag(z))/cmap->dy;
	int ix = (int)floor(wx); //Indices of the top left corner.
	int iy = (int)floor(wy);
	if (ix > cmap->nx-2) {//Clamp to existing indices.
		ix = cmap->nx-2;
	}
	else if (ix < 0) {
		ix = 0;
	}
	if (iy > cmap->ny-2) {
		iy = cmap->ny-2;
	}
	else if (iy < 0) {
		iy = 0;
	}
	wx -= ix; //Deduce the relative weights in [0, 1] of the four samples.
	wy -= iy;
	return (1-wx)*(1-wy)*get_value_cmap(cmap, ix, iy) + wx*(1-wy)*get_value_cmap(cmap, ix+1, iy) + (1-wx)*wy*get_value_cmap(cmap, ix, iy+1) + wx*wy*get_value_cmap(cmap, ix+1, iy+1);
}

/**
 * Computes the Laplacian of the given complex map "cmap" with the given kernel, and stores the result to the empty complex map of same dimensions "laplmap".
 * The pointer "laplmap" is supposed to be already allocated. Do not forget to call del_complexmap() to free the memory after use.
 */
void init_laplacemap(ComplexMap* cmap, ComplexMap* laplmap) {
	//if (cmapout->nx != cmap->nx && cmapout->ny != cmap->ny) {
	//	printf("[ERROR] The complex maps have not the same dimensions (%d,%d) != (%d,%d), aborting...\n", cmap->nx, cmap->ny, cmapout->nx, cmapout->ny); exit(1);}
	//if (cmapout->xmin != cmap->xmin && cmapout->xmax != cmap->xmax && cmapout->ymin != cmap->ymin && cmapout->ymax != cmap->ymax) {
	//	printf("[ERROR] The complex maps does not cover the same region of the complex plane, aborting...\n"); exit(1);}
	init_complexmap(laplmap, cmap->nx, cmap->ny, cmap->xmin, cmap->xmax, cmap->ymin, cmap->ymax, cmap->nseed, cmap->title);
	laplmap->realtime = cmap->realtime;
	*(laplmap->colrule) = *(cmap->colrule); //Copy pointed structures (content, not adresses).
	int ix, iy, cx, cy, l = 0;
	for (iy = 0; iy < cmap->ny; iy++) {
		for (ix = 0; ix < cmap->nx; ix++) {
			//Center of the second derivative in the horizontal x direction:
			cx = ix;
			if (cx == 0) {
				cx++;
			}
			else if (cx == cmap->nx-1) {
				cx--;
			}
			//Center of the second derivative in the vertical y direction:
			cy = iy;
			if (cy == 0) {
				cy++;
			}
			else if (cy == cmap->ny-1) {
				cy--;
			}
			laplmap->data[l] =
					(get_value_cmap(cmap, cx+1, iy) - 2*get_value_cmap(cmap, cx, iy) + get_value_cmap(cmap, cx-1, iy))/(cmap->dx*cmap->dx) 
					+ (get_value_cmap(cmap, ix, cy+1) - 2*get_value_cmap(cmap, ix, cy) + get_value_cmap(cmap, ix, cy-1))/(cmap->dy*cmap->dy);
			
			l++;
		}
	}
}

/**
 * Prints a short summary of the parameters used in the complex map for the standard output.
 */
void print_param_complexmap(ComplexMap* cmap, int verbose) {
	printf("[INFO] Complex map entitled '%s' of size %dx%d in region x=%g:%g y=%g:%g for %d seeds.\n", cmap->title, cmap->nx, cmap->ny, cmap->xmin, cmap->xmax, cmap->ymin, cmap->ymax, cmap->nseed);
	if (verbose) {
		printf("[INFO] Region aspect ratio is (ymax-ymin)/(xmax-xmin)=%g, pixmap anisotropy is dy/dx=%g. Computed in %gs.\n", (cmap->ymax-cmap->ymin)/(cmap->xmax-cmap->xmin), cmap->dy/cmap->dx, cmap->realtime);
	}
}

/**
 * Displays the first sxs elements of the given complex map to standard output.
 * The ordering is natural.
 * @note This function is intended to testing purpose only.
 */
void print_complexmap(ComplexMap* cmap, int s) {
	int sx = s;
	int sy = s;
	if (sx > cmap->nx) {
		sx = cmap->nx;
	}
	if (sy > cmap->ny) {
		sy = cmap->ny;
	}
	int ix, iy, l = 0;
	for (iy = 0; iy < sy; iy++) {
		for (ix = 0; ix < sx; ix++) {
			printf("\t%.3f%+.3fi", creal(cmap->data[l]), cimag(cmap->data[l]));
			l++;
		}
		printf("\n");
	}
}

/**
 * Saves the complex map data, in addition to other metadata, to the file titled "fname".
 */
void save_complexmap(ComplexMap* cmap, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[kplane]\ntitle=%s\n", cmap->title);  //Title of the complex map.
	fprintf(fp, "realtime=%g\n", cmap->realtime); //Total physical duration of the computation of the complex map (in seconds).
	fprintf(fp, "nseed=%d\n", cmap->nseed);
	fprintf(fp, "xsample=%d\n", cmap->nx);
	fprintf(fp, "ysample=%d\n", cmap->ny);
	fprintf(fp, "xrange=%g:%g\n", cmap->xmin, cmap->xmax);
	fprintf(fp, "yrange=%g:%g\n", cmap->ymin, cmap->ymax);
	save_colorrule(cmap->colrule, "color", fp);
	save_cdata((ComplexData*)cmap, fp);
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Plots the real part of the given complex map "cmap" to a PPM file (Portable PixMap).
 * This function then attempts to convert the PPM file to a PNG file if the appropriate command "pnmtopng" exists.
 * Note that the filename "fname" should not have any extension.
 */
void export_image_realpart(ComplexMap* cmap, int id, const char* fname) {
	size_t imgfnamelen = strlen(fname) + (size_t)floor(1+log10(1+abs(id))) + 2;
	char imgfname[imgfnamelen];
	sprintf(imgfname, "%s_%d", fname, id);
	find_min_max_cmap(cmap, 1); //Find minimum and maximum real parts and store them in the color rule (1=verbose).
	char ppmfile[strlen(imgfname)+6];
	sprintf(ppmfile, "%s.ppm", imgfname);
	FILE* outfp = fopen(ppmfile, "wb");  //'wb' for binary writing.
	fprintf(outfp, "P6\n%d %d\n255\n", cmap->nx, cmap->ny);  //Header with 'P6' for binary PPM file format.
	Pixel pix;
	int l;
	for (l = 0; l < cmap->ntot; l++) {//Loop on pixels to write.
		//set_color(creal(cmap->data[l]), minreal, maxreal, colorscheme, &pix); //Very old version
		//set_color(cmap->cols, creal(cmap->data[l]), minreal, maxreal, &pix); //Old version
		set_color(cmap->colrule, creal(cmap->data[l]), &pix); //This assumes "hmin" and "hmax" are already initialized (see find_min_max_cmap).
		fwrite(&pix, 1, sizeof(Pixel), outfp);
	}
	fclose(outfp);
	convert_ppm_to_png(imgfname);
}

/**
 * Appends a TikZ picture environment importing the image "fname.png" to a file "tikzfp".
 * Add the metadata in "cmap" to the plot, and a color legend.
 * Note that the filename "fname" should not have any extension.
 * @version In this version, the resonances of the square well approximation are drawn in the same axis as the complex map.
 */
void export_tikz_cmap(ComplexMap* cmap, Medium* med, int lmax, int id, const char* fname) {
	export_image_realpart(cmap, id, fname);
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	double txmin = cmap->xmin - cmap->dx/2; //True TikZ min/max ranges.
	double txmax = cmap->xmax + cmap->dx/2; //Indeed, the bounds of the complex map consider the center of each pixels.
	double tymin = cmap->ymin - cmap->dy/2;
	double tymax = cmap->ymax + cmap->dy/2;
	fprintf(tikzfp, "\\par%%\n\\begin{tikzpicture}[%%\n\tevery node/.append style={font=\\footnotesize},\n\t/pgfplots/width=250pt,\n\t/pgfplots/xlabel={$\\Re k~(\\varsigma^{-1})$},\n\t/pgfplots/ylabel={$\\Im k~(\\varsigma^{-1})$},\n\t/pgfplots/title style={align=center},\n\t/pgfplots/enlargelimits=false,\n\t/pgfplots/yticklabel style={rotate=90},\n\t/pgfplots/axis on top=true,\n\t/pgfplots/legend pos=outer north east,\n\t/pgfplots/legend columns=2,\n\tevery even column/.append style={column sep=5pt},\n\t/pgfplots/every node near coord/.style={font=\\tiny, inner sep=1.7pt},\n]%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\tname={cmap},\n\ttitle={%dD %s for $N=%d$ and %s\\\\ with %d seeds (map $%d\\times %d$, %.1f s)},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, cmap->nseed, cmap->nx, cmap->ny, cmap->realtime); //\txmin=%g, xmax=%g,\n\tymin=%g, ymax=%g,\n
	if (fabs(log(cmap->dy/cmap->dx)) < 0.2) {//If the pixel anisotropy is close enough to 1, then force the unit aspect ratio of the XY axes scale.
		fprintf(tikzfp, "\taxis equal image,\n");
	}
	else {//For large pixel anisotropy, ensure that the bitmap image is not distorted by the difference of scales of the XY axes.
		fprintf(tikzfp, "\tunit vector ratio*={1 %g},\n", cmap->dx/cmap->dy);
	}
	find_min_max_cmap(cmap, 0); //Find maximum value (0=quiet mode).
	char* colstr = tikz_colorbar_v2(cmap->colrule); //Colorbar v1=equally spaced ticks, v2=ticks are powers of ten.
	fprintf(tikzfp, "%s", colstr);
	free(colstr);
	fprintf(tikzfp, "]%%\n\\addplot[forget plot] graphics[xmin=%g, xmax=%g, ymin=%g, ymax=%g] {%s_%d.png};\n", txmin, txmax, tymin, tymax, fname, id);
	if (lmax >= 0) {//If lmax is positive, then also plots the resonances of the square well approximation using the medium "med" and the domain "dom".
		if (med->shape != ball) {//Warn the user if the medium is not ball shaped.
			printf("[WARN] The medium shape '%s' is not 'ball'. The resonances are not necessarily meaningful...\n", shape_to_string(med->shape));
		}
		printf("[INFO] Computing resonances of the effective square well approx...\n");
		Domain dom = {//Creates the domain, or create internal object of "cmap"?
			.xmin = cmap->xmin,
			.xmax = cmap->xmax,
			.ymin = cmap->ymin,
			.ymax = cmap->ymax,
		};
		char* str = resonance_plot_string(med, &dom, lmax); //Computes the resonances and gets the string.
		fprintf(tikzfp, "%s", str);
		free(str);
	}
	fprintf(tikzfp, "\\end{axis}\n\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * Appends a TikZ picture environment importing the image "fname.png" to a file "fname.tikz".
 * Add the metadata in "cmap" to the plot, and a color legend.
 * Note that the filename "fname" should not have any extension.
 * @version In this version, the resonances of the square well approximation are drawn in a separate axis.
 */
void export_tikz_cmap_v2(ComplexMap* cmap, Medium* med, int lmax, int id, const char* fname) {
	export_image_realpart(cmap, id, fname);
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	fprintf(tikzfp, "\\par%%\n\\begin{tikzpicture}[%%\n\tevery node/.append style={font=\\footnotesize},\n\t/pgfplots/width=250pt,\n\t/pgfplots/xmin=%g,\n\t/pgfplots/xmax=%g,\n\t/pgfplots/ymin=%g,\n\t/pgfplots/ymax=%g,\n\t/pgfplots/xlabel={$\\Re k~(\\varsigma^{-1})$},\n\t/pgfplots/ylabel={$\\Im k~(\\varsigma^{-1})$},\n\t/pgfplots/title style={align=center},\n\t/pgfplots/axis on top=true,\n\tevery even column/.append style={column sep=5pt},\n]%%\n\\begin{axis}[%%\n",
		cmap->xmin - cmap->dx/2, cmap->xmax + cmap->dx/2, cmap->ymin - cmap->dy/2, cmap->ymax + cmap->dy/2);
	//fprintf(tikzfp, "\ttitle={``%s''\\\\ with %d seeds (map $%d\\times %d$ computed in %g s)},\n\ttitle style={align=center},\n", cmap->title, cmap->nseed, cmap->nx, cmap->ny, cmap->realtime);  //OLD version
	fprintf(tikzfp, "\tname={cmap},\n\ttitle={%dD %s for $N=%d$ and %s\\\\ with %d seeds (map $%d\\times %d$, %.1f s)},\n",
		med->d, shape_to_string(med->shape), med->n, med->scmodel->param_latex, cmap->nseed, cmap->nx, cmap->ny, cmap->realtime);
	if (lmax >= 0) {//If lmax positive, then hide the ticks and labels to save place for subplot below.
		fprintf(tikzfp, "\txlabel=\\empty,\n\txticklabels=\\empty,\n");
	}
	fprintf(tikzfp, "\tyticklabel style={rotate=90},\n\taxis equal image,\n");
	find_min_max_cmap(cmap, 0); //Find maximum value (0=quiet mode).
	char* colstr = tikz_colorbar_v2(cmap->colrule);
	fprintf(tikzfp, "%s", colstr);
	free(colstr);
	fprintf(tikzfp, "]%%\n\\addplot graphics[xmin=\\pgfkeysvalueof{/pgfplots/xmin}, xmax=\\pgfkeysvalueof{/pgfplots/xmax}, ymin=\\pgfkeysvalueof{/pgfplots/ymin}, ymax=\\pgfkeysvalueof{/pgfplots/ymax}] {%s_%d.png};\n\\end{axis}\n", fname, id);
	if (lmax >= 0) {//If lmax is positive, then also plots the resonances of the square well approximation using the medium "med" and the domain "dom".
		printf("[INFO] Computing resonances of the effective square well approx...\n");
		Domain dom = {//Creates the domain, or create internal object of "cmap"?
			.xmin = cmap->xmin,
			.xmax = cmap->xmax,
			.ymin = cmap->ymin,
			.ymax = cmap->ymax,
		};
		char* str = resonance_plot_string(med, &dom, lmax);  //Computes the resonances and gets the string.
		fprintf(tikzfp, "\\begin{axis}[%%\n\tname={sqwell-reson-approx},\n\tat={(cmap.south west)},\n\tanchor={north west},\n\tyshift=-10pt,\n\tlegend pos={outer north east},\n\tlegend columns=2,\n\textra x ticks={0},\n\textra y ticks={0},\n\textra tick style={grid=major},\n\textra x tick labels=\\empty,\n\textra y tick labels=\\empty,\n\tyticklabel style={rotate=90},\n\taxis equal image,\n]%%\n");
		fprintf(tikzfp, "%s", str);
		fprintf(tikzfp, "\\end{axis}\n");
		free(str);
	}
	fprintf(tikzfp, "\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * Appends a plot of the vertical cut of the given complex map "cmap" at the real part "x" to the tikzpicture file "fname.tikz".
 * The abscissa is -Im(k), and thus extends from -ymax to -ymin. The samples are obtained through bilinear interpolation.
 */
void export_tikz_cut(ComplexMap* cmap, double x, const char* fname) {
	char tikzfile[strlen(fname)+6];
	sprintf(tikzfile, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfile, "a");
	fprintf(tikzfp, "\\par\\begin{tikzpicture}\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\txmin=%g, xmax=%g,\n", -cmap->ymax, -cmap->ymin);
	fprintf(tikzfp, "\txlabel={$-\\Im k~(\\varsigma^{-1})$},\n\tylabel={%s},\n", cmap->title);
	fprintf(tikzfp, "\ttitle={Vertical cut at $\\Re k=%g$},\n", x);
	fprintf(tikzfp, "\twidth=0.4\\textwidth,\n");
	fprintf(tikzfp, "]\n\\addplot[red, thick, mark=none] coordinates {%%\n");
	double y;
	dcomplex z, fz;
	int i, n = 150; //Number of samples.
	for (i = 0; i < n; i++) {//Loop on samples along the vertical cut.
		y = cmap->ymax + (cmap->ymin-cmap->ymax)*i/(n-1);
		z = x + y*I;
		fz = interpolate_bilinear(cmap, z);
		fprintf(tikzfp, "(%g, %g) ", -y, creal(fz)); //TikZ coordinates format.
	}
	fprintf(tikzfp, "\n};\n\\end{axis}\n\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * Parses the complex map "cmap" with the arguments of "args".
 * The "cmap" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_complexmap(ComplexMap* cmap, int narg, char** args) {
	int nargmin = 4; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(1);
	}
	double xmin, xmax, ymin, ymax;
	int nx = 0;  //Default value.
	int ny = 0;  //Default value.
	int nseed;
	char *title, *value1, *value2;
	value1 = get_value(narg, args, "xrange");
	if (sscanf(value1, "%lg:%lg", &xmin, &xmax) != 2) {
		printf("[ERROR] Invalid real range '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "yrange");
	if (sscanf(value1, "%lg:%lg", &ymin, &ymax) != 2) {
		printf("[ERROR] Invalid imaginary range '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "xsample");
	value2 = get_value(narg, args, "ysample");
	int match1 = sscanf(value1, "%d", &nx); //If done in a condition statement, the optimizer will skip the second sscanf().
	int match2 = sscanf(value2, "%d", &ny);
	if (match1 != 1 && match2 != 1) {//Just need one number of sample, the other one can be inferred.
		printf("[ERROR] Invalid sample numbers nx='%s', ny='%s', aborting...\n", value1, value2);
		exit(1);
	}
	value1 = get_value(narg, args, "nseed");
	if (sscanf(value1, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value1);
		exit(1);
	}
	title = get_value(narg, args, "title");
	init_complexmap(cmap, nx, ny, xmin, xmax, ymin, ymax, nseed, title); //Call initialization function.
	sscanf(get_value(narg, args, "realtime"), "%lg", &(cmap->realtime)); //Set realtime after initialization.
	value1 = get_value(narg, args, "color"); //Now parsing the color scheme.
	parse_colorrule(cmap->colrule, value1);
	value1 = get_value(narg, args, "data");  //Now parsing the list of data.
	parse_cdata((ComplexData*)cmap, value1);
}
