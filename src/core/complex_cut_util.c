/****
 * @date Created on 2020-11-30 at 14:40:43 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities for manipulating complex cuts.
 ***/
#include <stdlib.h>             /* Standard Library for Memory Allocation */
#include <stdio.h>              /* Standard Library for Input and Output */
#include <string.h>             /* Standard Library for String Manipulation */
#include "common_util.h"        /* Import General Purpose Functions such as the Progress Bar */
#include "complex_cut_type.h"   /* Import the Complex Cut Type */
#include "complex_data_util.h"  /* Import the Complex Data Parsing functions */
#include "medium_type.h"        /* Import the Medium Type */
/**
 * Define some constant macros:
 */
#define MEPS     1.0e-16        /* Approximate machine epsilon of double precision. */
#define LC0      -3.0           /* Discrete Laplacian coefficient for the central component. */
#define LC1       0.5           /* Discrete Laplacian coefficient for the nearest component. */
#define LC2       0.25          /* Discrete Laplacian coefficient for the diagonal component. */

/**
 * Returns the position in the z-plane corresponding to the given index "l" in the range {0, ..., 3*ns+5}.
 * The samples are assumed to be stored in triplet-major format, running in the transverse direction from left to right,
 * then in the longitundinal direction from "zmin" to "zmax".
 */
dcomplex get_arg_ccut(ComplexCut* ccut, int l) {
	int il = l/3;  //Longitudinal index coordinate.
	int it = l%3;  //Transverse index coordinate.
	return ((ccut->ns-il)*ccut->zmin + (il-1)*ccut->zmax)/(ccut->ns-1) + I*ccut->dz*(it-1);
}

/**
 * Returns the element of the complex cut at the given index "i" in the range {0, ..., ns-1}.
 * The origin is at "zmin". Do not check for bounds.
 */
dcomplex get_value_ccut(ComplexCut* ccut, int i) {
	return ccut->data[3*i+4];
}

/**
 * Returns the Laplacian of the complex cut at the given index "i" in the range {0, ..., ns-1}.
 * Use a nine-point stencil. Do not check for bounds.
 */
dcomplex get_laplacian_ccut(ComplexCut* ccut, int i) {
	return (LC2*ccut->data[3*i] + LC1*ccut->data[3*i+1] + LC2*ccut->data[3*i+2] + LC1*ccut->data[3*i+3]
		+ LC0*ccut->data[3*i+4] + LC1*ccut->data[3*i+5] + LC2*ccut->data[3*i+6] + LC1*ccut->data[3*i+7] + LC2*ccut->data[3*i+8])/(ccut->h*ccut->h);
}

/**
 * Initialization function of an empty complex cut with the given parameters.
 * The content of the complex cut must be deleted with the corresponding del_*() function.
 */
void init_complexcut(ComplexCut* ccut, int ns, dcomplex zmin, dcomplex zmax, int nseed, const char* title) {
	if (ns > 10000000) {//Beyond 10 MegaPixels is a bit too large...
		printf("[ERROR] The number of samples is too large (ns=%d), aborting...\n", ns);
		exit(1);
	}
	if (ns < 2) {
		printf("[ERROR] The number of samples is too small (ns=%d), aborting...\n", ns);
		exit(1);
	}
	dcomplex dz = (zmax-zmin)/(ns-1);
	double h = cabs(dz);
	if (h < MEPS*cabs(zmin) || h < MEPS*cabs(zmax)) {
		printf("[ERROR] The end point is too close to the start point (step |dz|=%g), aborting...\n", h);
		exit(1);
	}
	ccut->ns = ns;
	ccut->ntot = 3*(ns+2);
	ccut->zmin = zmin;
	ccut->zmax = zmax;
	ccut->dz = dz;
	ccut->h = h;
	ccut->nseed = nseed;
	strncpy(ccut->title, title, 79); //Maximum size of the title.
	ccut->realtime = 0.0; //Computation time is zero, since empty.
	ccut->data = (dcomplex*)calloc(ccut->ntot, sizeof(dcomplex));
	ccut->get_arg = (dcomplex(*)(void*, int))get_arg_ccut;
}

/**
 * Deletes the content (allocated pointers) within the given complex cut, not the pointer itself.
 * Do not forget to invoke after use of "init_*" functions.
 */
void del_complexcut(ComplexCut* ccut) {
	free(ccut->data);
}

/**
 * Prints a short summary of the parameters used in the complex cut for the standard output.
 */
void print_param_complexcut(ComplexCut* ccut) {
	printf("[INFO] Complex cut entitled '%s' with %d samples (%dx%d) and endpoints zmin=%g%+gi and zmax=%g%+gi for %d seeds.\n",
		ccut->title, ccut->ns, ccut->ns+2, 3, creal(ccut->zmin), cimag(ccut->zmin), creal(ccut->zmax), cimag(ccut->zmax), ccut->nseed);
}

/**
 * Saves the complex cut data, in addition to other metadata, to the file titled "fname".
 */
void save_complexcut(ComplexCut* ccut, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[kcut]\ntitle=%s\n", ccut->title);  //Title of the complex cut.
	fprintf(fp, "realtime=%g\n", ccut->realtime); //Total physical duration of the computation of the complex cut (in seconds).
	fprintf(fp, "nseed=%d\n", ccut->nseed);
	fprintf(fp, "sample=%d\n", ccut->ns);
	fprintf(fp, "kmin=%.12g%+.12gi\n", creal(ccut->zmin), cimag(ccut->zmin));
	fprintf(fp, "kmax=%.12g%+.12gi\n", creal(ccut->zmax), cimag(ccut->zmax));
	save_cdata((ComplexData*)ccut, fp);
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Appends a plot of the complex cut to the tikzpicture file "fname.tikz".
 * If "lflag" is 1, then computes the Laplacian of the complex cut data.
 */
void export_tikz_ccut(ComplexCut* ccut, Medium* med, int lflag, const char* fname) {
	char tikzfile[strlen(fname)+6];
	sprintf(tikzfile, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfile, "a");
	fprintf(tikzfp, "\\par\\begin{tikzpicture}\n\\begin{axis}[%%\n");
	//fprintf(tikzfp, "\txmin=0, xmax=%g,\n", cabs(ccut->zmax-ccut->zmin));  //Length travelled inthe k-plane (in units of "k").
	fprintf(tikzfp, "\txlabel={Coordinate ($\\varsigma^{-1}$)},\n\tylabel={%s},\n", ccut->title);
	fprintf(tikzfp, "\ttitle={Cut from $%g%+g\\,{\\rm i}$ to $%g%+g\\,{\\rm i}$\\\\ with %d seeds (computed in %g s)},\n\ttitle style={align=center},\n", creal(ccut->zmin), cimag(ccut->zmin), creal(ccut->zmax), cimag(ccut->zmax), ccut->nseed, ccut->realtime);
	if (lflag) {//Enable bilogarithmic axes only for resonance density.
		fprintf(tikzfp, "\tymode=log,\n"); //xmode=log, 
		fprintf(tikzfp, "\tlegend entries={Computation,$(d+3)/4k_{\\rm i}^2$},\n\tlegend pos=south west,\n");
	}
	fprintf(tikzfp, "\twidth=0.6\\textwidth,\n");
	fprintf(tikzfp, "]\n\\addplot[black, thick, mark=none] coordinates {%%\n\t");
	double x;
	dcomplex fz;
	int i;
	for (i = 0; i < ccut->ns; i++) {//Loop on the effective samples (central samples).
		x = i*ccut->h;  //Longitudinal coordinate (in units of "k").
		if (lflag) {
			fz = get_laplacian_ccut(ccut, i);
		}
		else {
			fz = get_value_ccut(ccut, i);
		}
		fprintf(tikzfp, "(%g, %g) ", x, creal(fz)); //TikZ coordinates format.
	}
	fprintf(tikzfp, "\n};\n\\addplot[red, thick, mark=none] coordinates {%%\n\t");
	if (lflag) {//Plots our best upper bound for the ball medium: (d+3)/(4*imk^2). Also an approx for other shapes.
		for (i = 0; i < ccut->ns; i++) {
			x = cimag(ccut->zmin) + cimag(ccut->zmax-ccut->zmin)*i/(ccut->ns-1);  //Longitudinal coordinate (in units of "k").
			if (MEPS < x || x < -MEPS) {
				fprintf(tikzfp, "(%g, %g) ", i*ccut->h, (med->d+3)/(4*x*x)); //TikZ coordinates format.
			}
		}
	}
	fprintf(tikzfp, "\n};\n\\end{axis}\n\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * Parses the complex cut "ccut" with the arguments of "args".
 * The "ccut" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_complexcut(ComplexCut* ccut, int narg, char** args) {
	int nargmin = 4; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(1);
	}
	int ns;
	double re, im;
	dcomplex zmin, zmax;
	int nseed;
	char *title, *value1;
	value1 = get_value(narg, args, "sample");
	if (sscanf(value1, "%d", &ns) != 1) {
		printf("[ERROR] Invalid number of samples '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "kmin");
	if (sscanf(value1, "%lg%lgi", &re, &im) != 2) {
		printf("[ERROR] Invalid endpoint '%s', aborting...\n", value1);
		exit(1);
	}
	zmin = re + im*I;
	value1 = get_value(narg, args, "kmax");
	if (sscanf(value1, "%lg%lgi", &re, &im) != 2) {
		printf("[ERROR] Invalid endpoint '%s', aborting...\n", value1);
		exit(1);
	}
	zmax = re + im*I;
	value1 = get_value(narg, args, "nseed");
	if (sscanf(value1, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value1);
		exit(1);
	}
	title = get_value(narg, args, "title");
	init_complexcut(ccut, ns, zmin, zmax, nseed, title); //Call initialization function.
	sscanf(get_value(narg, args, "realtime"), "%lg", &(ccut->realtime)); //Set realtime after initialization.
	value1 = get_value(narg, args, "data");  //Now parsing the list of data.
	//parse_complexcut_data(ccut, value1);
	//parse_complex_data(ccut->ntot, ccut->data, value1);
	parse_cdata((ComplexData*)ccut, value1);
}

