/****
 * @date Created on 2021-02-18 at 13:59:41 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing utilities to generate TikZ/PGFPlots figures from data.
 ***/
#include <stdlib.h>        /* Standard Library for Memory Allocation */
#include <stdio.h>         /* Standard Library for Input and Output */
#include <string.h>        /* Standard Library for String Manipulation */
#include "curve_type.h"    /* Imports the Curve Type */

/**
 * Allocates and initializes the curve "crv".
 */
void init_curve(Curve* crv, int nx, const char* opt) {
	crv->nx = nx;
	crv->x = (double*)calloc(nx, sizeof(double));
	crv->y = (double*)calloc(nx, sizeof(double));
	strncpy(crv->opt, opt, 79);
}

/**
 * Removes the content of the curve "crv".
 */
void del_curve(Curve* crv) {
	free(crv->x);
	free(crv->y);
}

/**
 * Writes in the file "fp" a TikZ "axis" environment to plot the given curves stored in "crv".
 * The number of curves to plot is "nc" (>=1). The file is not closed.
 * The string "axisopt" contains additional options for the "axis" environment.
 * The string "precmd" appends possible additional TikZ commands before the addplot command in the axis environment (must end with line return).
 */
void tikz_plot(FILE* fp, Curve* crv, int nc, char* axisopt, char* precmd) {
	fprintf(fp, "\\begin{axis}[%%\n\t%s\n]\n%s", axisopt, precmd);
	int i, j;
	for (j = 0; j < nc; j++) {//Loop on the curves to plot.
		fprintf(fp, "\\addplot[%s] coordinates {%%\n\t", crv[j].opt);
		for (i = 0; i < crv[j].nx; i++) {//Loop on the data points.
			fprintf(fp, "(%.15lg, %.15lg) ", crv[j].x[i], crv[j].y[i]);
		}
		fprintf(fp, "\n};\n");
	}
	fprintf(fp, "\\end{axis}\n");
}

/**
 * Writes in the file "fp" a TikZ "axis" environment to plot the given points represented by the abscissa "xdata" of length "nx", and the ordinates "ydata" or length "nx*ny".
 * The "xdata" is stored in reading order (x_1, x_2, x_3, ...), and "ydata" in the same order (y_1, y_2, ..., y_nx)_plot1, (y_1, y_2, ..., y_nx)_plot2, etc.
 * The number of abscissa points is "nx" and the number of curves to plot is "ny" (>=1). The file is not closed.
 * The string "axisopt" contains additional options for the "axis" environment,
 * and the array of strings "plotopt" of length "ny" contains the options for each curve to be drawn.
 * @deprecated Old version which does not allow to plot curves with different numbers of points.
 */
void tikz_plot_old(FILE* fp, double* xdata, int nx, double* ydata, int ny, char* axisopt, char** plotopt) {
	fprintf(fp, "\\begin{axis}[%%\n");
	fprintf(fp, "\t%s,\n", axisopt);
	fprintf(fp, "\tyticklabel style={rotate=90},\n]\n");
	int i, j;
	for (j = 0; j < ny; j++) {//Loop on the curves to plot.
		fprintf(fp, "\\addplot[%s] coordinates {%%\n", plotopt[j]);
		for (i = 0; i < nx; i++) {//Loop on the data points.
			fprintf(fp, "(%.15lg, %.15lg) ", xdata[i], ydata[i + j*nx]);
		}
		fprintf(fp, "\n};\n");
	}
	fprintf(fp, "\\end{axis}\n");
}

/**
 * Same as tikz_plot() but using the explicit interval of abscissas [xmin, xmax].
 * @deprecated Unused because simpler to directly construct both "xdata" and "ydata" at the same time when needed.
 */
void tikz_plot_range(FILE* fp, double xmin, double xmax, int nx, double* ydata, int ny, char* axisopt, char** plotopt) {
	double dx = (xmax - xmin)/(nx - 1);
	double* xdata = (double*)calloc(nx, sizeof(double));
	int i;
	for (i = 0; i < nx; i++) {
		xdata[i] = xmin + i*dx;
	}
	tikz_plot_old(fp, xdata, nx, ydata, ny, axisopt, plotopt);
	free(xdata);
}
