/****
 * @date Created on 2021-02-01 at 09:27:32 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to solve the pole of the classical transport problem in a uniform spherical medium in dimension "d".
 * Use a finite volume method to discretize the integral transport equation in the d-ball,
 * and a Rayleigh iteration to find-root the pole in "beta" (sp^-1) and the transport eigenmode.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "common_util.h"
#define umi(i, j) (i + (j)*(j + 1)/2)  //Define the absolute index of packed upper triangle matrix in terms of the std indices (i,j). We must have i <= j (!).

//LAPACK's routine prototype for real packed symmetric matrices:
void dspev_(char* jobz, char* uplo, int* n, double* ap, double* w, double* z, int* ldz, double* work, int* info);

typedef struct BallGrid_s {
	int d;         //Number of spatial dimensions.
	int n;         //Number of atoms in the quantum case.
	double sigma;  //Total cross section of single atoms in the quantum case (units of sp^(d-1)). Since the density is unity, also the inverse of the mean free path.
	double radius; //Radius of the ball medium from the quantum case, ensuring unit density (units of sp).
	int s;         //Number of grid points along a diameter of the ball, so that s^d is the maximmum number of grid points (cubic upper bound).
	double h;      //Distance between two consecutive grid points (units of sp). Step of the cubic lattice.
	double v;      //Mean volume of one of the cubic cells, equal to h^d.
	int np;        //Number of grid points in "pos". The array "pos" generally contains more points (s^d), but the last coordinates are left empty.
	double* pos;   //Positions of the grid points, stored in natural reading order (x1,y1,z1) (x2,y2,z2), etc. totalling d*s^d coordinates.
} BallGrid;

/**
 * Initializes the given grid points array "bg" with at most s^d points inside the ball of given radius.
 * The integer "s" represents the number of grid points along the diameter of the ball.
 * The points outside of the ball are not stored. Returns the number of grid points in the ball.
 */
void init_ball_grid(BallGrid* bg, int d, int n, double sigma, int s) {
	bg->d = d;
	bg->n = n;
	bg->sigma = sigma;
	bg->s = s;
	bg->radius = pow(n/uball_volume(d), 1./d);  //Radius of the d-ball to ensure unit density with "n" atoms.
	bg->h = (2.*bg->radius)/s;  //Distance between consecutive grid points.
	bg->v = pow(bg->h, d);      //Typical volume of a cubic cell.
	int pmax = ipow(s, d);      //Maximum number of grid points in a cube of side "s".
	bg->pos = (double*)calloc(d*pmax, sizeof(double)); //Allocate memory for the grid points (totalling d*s^d coordinates).
	int p, val, c, l, i = 0; //i=Index of the currently processed grid points.
	double norm2, r2 = bg->radius*bg->radius;
	for (p = 0; p < pmax; p++) {//Loop over the grid points of the cubic lattice.
		val = p;
		l = i*d; //Absolute index in the array "pos".
		norm2 = 0.0;
		for (c = 0; c < d; c++) {//Change of base method ("val" in base "s").
			bg->pos[l+c] = bg->h*(val%s - (s-1.)/2);
			val /= s;
			norm2 += bg->pos[l+c]*bg->pos[l+c];  //Compute the distance from the center.
		}
		if (norm2 < r2) {//Skip to the next point only if the current grid point lies in the ball, otherwise the point is overwritten.
			i++;
		}
	}
	bg->np = i;  //Stores the actual number of grid points in the ball.
}

/**
 * Deletes the content (allocated pointers) in the given ball grid, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_ball_grid(BallGrid* bg) {
	free(bg->pos);
}

/**
 * Prints the given ball grid "bg" for checking purpose only.
 */
void print_ball_grid(BallGrid* bg, int verbose) {
	printf("[INFO] Ball Grid with dim = %d, s = %d, radius = %g sp, np = %d, h = %g sp.\n", bg->d, bg->s, bg->radius, bg->np, bg->h);
	if (verbose) {
		int i, c;
		for (i = 0; i < bg->np; i++) {
			printf("#%d = ", i);
			for (c = 0; c < bg->d; c++) {
				printf("\t%.4g", bg->pos[i*bg->d+c]);
			}
			printf("\n");
		}
	}
	fflush(stdout);
}

/**
 * Computes the Euclidean distance between the points of index "i" and "j" in the given list of coordinates "pos".
 * This function does not perform bounds checking.
 */
double cl_distance(double* pos, int i, int j, int d) {
	double dx, dist2 = 0.0;
	int c;
	for (c = 0; c < d; c++) {
		dx = pos[i*d + c] - pos[j*d + c];
		dist2 += dx*dx;
	}
	return sqrt(dist2);
}

/**
 * Estimates the diagonal entry of the system matrix using a heuristic approximation.
 * This function only works for the specific dimensions 1D, 2D, 3D, 4D, according to the limitation of the array containing the effective radii.
 */
static double nuls[] = {0.715, 0.816, 0.965};  //Unitless effective radius of cubic cell in 2D, 3D, 4D. Heuristic formula not reliable in 4D.
double matrix_diag_approx(BallGrid* bg, double beta) {
	//diag = (beta + bg->sigma*exp(-(beta+bg->sigma)*bg->h/2.5))/(beta + bg->sigma);  //OLD: Heuristics to mimic the integral diagonal term Integral[K(p1-p2) , p1, p2] in the same cell.
	double ita;  //Approx of: (1/V(C))*Integral[exp(-gamma*|x1-x2|)/Surf(d)*|x1-x2|^(d-1), x1 in C, x2 in C], where C=[0,h]^d denotes a d-cube of side "h". Units of "sp" (inter-atomic spacing).
	double gamma = beta + bg->sigma;
	if (bg->d == 1) {//Exact result in 1D.
		double g = gamma*bg->h;
		ita = (g - 1 + exp(-g))/(gamma*g);
	}
	else {//Heuristic estimate in higher dimensions.
		double nu = nuls[bg->d-2];  //Unitless effective radius.
		ita = (1 - exp(-gamma*nu*bg->h))/(uball_volume(bg->d)*pow(nu, bg->d)*gamma);
	}
	return 1. - bg->sigma*ita;  //Diagonal entry approximation.
}

/**
 * Builds the classical multiple scattering matrix "ap" in the discretized ball-shaped medium "bg" with time parameter "beta".
 * Use packed column-major format to store only the upper triangular part in "ap". So, "ap" has length np*(np+1)/2.
 */
void build_cl_matrix(BallGrid* bg, double beta, double* ap) {
	double dist, diag = matrix_diag_approx(bg, beta);  //Approximates the diagonal entry of the matrix more precisely.
	//printf("[INFO] Setting diag = %g \n", diag); fflush(stdout);
	double fac = bg->sigma*bg->v/uball_surface(bg->d);  //Pre-computes the common prefactor of the off-diagonal elements.
	int i, j;
	for (j = 0; j < bg->np; j++) { //Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) { //Loop over the rows in the upper triangle (matrix is symmetric).
			dist = cl_distance(bg->pos, i, j, bg->d);
			ap[umi(i, j)] = -fac*exp(-(beta+bg->sigma)*dist)/pow(dist, bg->d-1);  //Off-diagonal term.
		}
		ap[umi(j, j)] = diag; //Diagonal term.
	}
}

/**
 * Prints to stdout the given real symmetric "n" x "n" matrix stored as "ap" in packed column-major format (upper triangle).
 * The array "ap" must contain at least n*(n+1)/2 elements.
 */
void print_packed_matrix(double* ap, int n, const char* matrixname) {
	int i, j;
	printf("[INFO] %s = \n", matrixname);
	for (i = 0; i < n; i++) {//Loop on rows.
		for (j = 0; j < n; j++) {//Loop on columns.
			if (j < i) {//Lower triangle.
				printf("\t%8.3g", ap[umi(j, i)]);
			}
			else {//Upper triangle.
				printf("\t%8.3g", ap[umi(i, j)]);
			}
		}
		printf("\n");
	}
	fflush(stdout);
}

/**
 * Function intended to print the matrix A for the given parameters.
 * d=Dimension, n=Number of atoms (quantum case), sigma=Total single-atom cross section also inverse of mean free path (quantum case),
 * beta=Time dual variable (sp^-1), s=Number of grid points along the ball diameter.
 */
void test_matrix(int d, int n, double sigma, int s, double beta) {
	BallGrid* bg = (BallGrid*)calloc(1, sizeof(BallGrid));
	init_ball_grid(bg, d, n, sigma, s);
	print_ball_grid(bg, 1);  //Prints a status for checking purpose.
	fflush(stdout);
	double* ap = (double*)calloc(bg->np*(bg->np + 1)/2, sizeof(double));
	printf("[INFO] Successfully allocating matrix with size %dx%d (%d elements)...\n", bg->np, bg->np, bg->np*(bg->np + 1)/2); fflush(stdout);
	build_cl_matrix(bg, beta, ap);
	print_packed_matrix(ap, bg->np, "System matrix A");
	//printf("[INFO] Minimum eigenvalue = %.15g \n", min_eigval(ap, bg->np));
	free(ap);
	del_ball_grid(bg);
	free(bg);
}

/**
 * Find the minimum element in absolute value of the given list "ls" of length "n".
 */
double find_min(double* ls, int n) {
	if (n == 0) {
		printf("[WARN] The list is empty, no minimum was found.\n");
		return 0.0;
	}
	int i;
	double minval = ls[0];
	for (i = 1; i < n; i++) {
		if (fabs(ls[i]) < fabs(minval)) {
			minval = ls[i];
		}
	}
	return minval;
}

/**
 * Computes and returns the eigenvalue of the given real symmetric "n" x "n" matrix "ap" (packed format, upper triangle) which is closest to zero.
 * Uses LAPACK's dspev() routine. The matrix "ap" is modified upon execution.
 */
double min_eigval(double* ap, int n) {
	int info;
	double w[n], work[3*n];  //Eigenvalues, and buffer space for LAPACK.
	
	dspev_("Novector", "Upper", &n, ap, w, NULL, &n, work, &info);
	
	if (info < 0) {
		printf("[LAPACK] Argument #%d has an illegal value.\n", -info);
	}
	else if (info > 0) {
		printf("[LAPACK] The eigenvalue algorithm failed to converge.\n");
	}
	return find_min(w, n);
}

typedef struct MinEigvalData {
	int ns;       //Number of different ball grids, typically differing by their number of grid points.
	BallGrid* bg; //List of initialized ball lattices used for the computations.
	double bmin;  //Minimum value of beta. First value scanned in line search loops.
	double bmax;  //Maximum value of beta. Last value scanned in line search loops.
	int nb;       //Total number of beta samples. Same for the "ns" different simulations.
	double* data; //Data of length nb*(ns+1) containing the minimum eigenvalues, stored in reading order (b mev1 mev2 ...) (b mev1 mev2 ...) ...
	// where "mev1" "mev2" corresponds to the different ball grids (typically different number of lattice grid points).
} minev;

/**
 * Computes the minimum eigenvalues of the system matrix as a function of "beta", i.e., min_eigval(beta), for the i^th simulation stored in "mv".
 * The index must be contained in [0, mv->ns-1].
 */
void run_minev_sim(minev* mv, int i) {
	BallGrid* b = mv->bg + i;
	double* ap = (double*)calloc(b->np*(b->np + 1)/2, sizeof(double));
	double beta;
	int j;
	for (j = 0; j < mv->nb; j++) {//Loop on the couples (beta, min_eigval(beta)) to compute.
		beta = mv->data[j*(mv->ns+1)];  //Get the "beta" value associated with the current sample.
		build_cl_matrix(b, beta, ap);
		mv->data[i + 1 + j*(mv->ns+1)] = min_eigval(ap, b->np);
	}
	free(ap);
}

/**
 * Runs all the simulations planned in the structure "mv".
 */
void run_minev(minev* mv) {
	double db = (mv->bmax - mv->bmin)/(mv->nb - 1);
	int i, j;
	for (j = 0; j < mv->nb; j++) {//Loop on the abscissas of the "beta" samples.
		mv->data[j*(mv->ns+1)] = mv->bmin + j*db;  //Values of "beta" from "bmin" to "bmax".
	}
	for (i = 0; i < mv->ns; i++) {//Loop on the different simulations.
		printf("[EXEC] Running %d/%d, with %d grid points...\n", i, mv->ns, (mv->bg+i)->np);
		run_minev_sim(mv, i);
	}
}

/**
 * Fills the given array "minevls" of length 2*nb with the couples (beta, min_eigval(beta)), i.e., the minimum eigenvalue from the system matrix A(beta).
 * The values of "beta" are taken from "bmin" to "bmax" with step (bmax-bmin)/(nb-1). The lattice "bg" is assumed to be already initialized.
 * @deprecated Replaced by run_minev()
 */
void min_eigval_array(BallGrid* bg, double bmin, double bmax, int nb, double* minevls) {
	double beta, db = (bmax-bmin)/(nb-1);
	double* ap = (double*)calloc(bg->np*(bg->np + 1)/2, sizeof(double));
	int i;
	for (i = 0; i < nb; i++) {//Loop on the couples (beta, min_eigval(beta)) to compute.
		beta = bmin + i*db;
		build_cl_matrix(bg, beta, ap);
		minevls[2*i] = beta;
		minevls[2*i+1] = min_eigval(ap, bg->np);
	}
	free(ap);
}

/**
 * Writes a TikZ file to plot the given points "xydata" of length np*(nc+1) in the file "fname" (without ".tikz" extension).
 * The data is stored in reading order (x, y1, y2, ..., ync), (x, y1, y2, ..., ync), ... for each sample of given "x".
 * The number of abscissa points is "np" and the number of curves to plot is "nc" (>=1).
 * @deprecated Replaced by tikz_plot_minev()
 */
void tikz_plot(double* xydata, int np, int nc, const char* title, const char* xlabel, const char* ylabel, const char* fname) {
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	prepare_file(tikzfname, "%Computed");
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	fprintf(tikzfp, "\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n");
	//fprintf(tikzfp, "\txmin=%g, xmax=%g,\n", xydata[0], xydata[(np-1)*(nc+1)]);
	fprintf(tikzfp, "\txlabel={%s},\n\tylabel={%s},\n", xlabel, ylabel);
	fprintf(tikzfp, "\ttitle={%s},\n\ttitle style={align=center},\n", title);
	fprintf(tikzfp, "\textra x ticks={%g, %g}, extra y ticks={0},\n\textra tick style={grid=major}, extra x tick labels={},\n", xydata[0], xydata[(np-1)*(nc+1)]);
	fprintf(tikzfp, "\twidth=0.7\\textwidth,\n\tyticklabel style={rotate=90},\n]\n");
	int i, j;
	for (i = 1; i <= nc; i++) {//Loop on the curves to plot.
		fprintf(tikzfp, "\\addplot[mark=*, mark size=0.5] coordinates {\n");
		for (j = 0; j < np; j++) {//Loop on the data points.
			fprintf(tikzfp, "(%g, %g) ", xydata[j*(nc+1)], xydata[i + j*(nc+1)]);
		}
		fprintf(tikzfp, "\n};\n");
	}
	fprintf(tikzfp, "\\end{axis}\n\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * TODO: Add legend to the plot.
 */
static const char* colors[] = {"red", "orange!80!black", "green!50!black", "cyan!80!blue", "blue"};
static const int ncolors = 5;
void tikz_plot_minev(minev* mv, const char* fname) {
	int i, j;
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	prepare_file(tikzfname, "%Computed");
	FILE* tikzfp = fopen(tikzfname, "a"); //Append to possibly existing TikZ file.
	fprintf(tikzfp, "\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n");
	fprintf(tikzfp, "\txlabel={$\\beta$~($\\varsigma^{-1}$)},\n\tylabel={Min. eigenvalue},\n");
	fprintf(tikzfp, "\ttitle={Minimum eigenvalue of $A(\\beta)$},\n\ttitle style={align=center},\n");
	fprintf(tikzfp, "\tlegend entries={");
	for (i = 0; i < mv->ns; i++) {//Loop on the curves to plot.
		fprintf(tikzfp, "$s=%d$,", (mv->bg+i)->s);
	}
	fprintf(tikzfp, "}, legend pos={south east},\n");
	fprintf(tikzfp, "\textra x ticks={%g, %g}, extra y ticks={0},\n\textra tick style={grid=major}, extra x tick labels={},\n", mv->data[0], mv->data[(mv->nb-1)*(mv->ns+1)]);
	fprintf(tikzfp, "\twidth=0.7\\textwidth,\n\tyticklabel style={rotate=90},\n]\n");
	for (i = 0; i < mv->ns; i++) {//Loop on the curves to plot.
		fprintf(tikzfp, "\\addplot[%s, mark=*, mark size=0.5] coordinates {\n", colors[i%ncolors]);
		for (j = 0; j < mv->nb; j++) {//Loop on the data points.
			fprintf(tikzfp, "(%g, %g) ", mv->data[j*(mv->ns+1)], mv->data[i + 1 + j*(mv->ns+1)]);
		}
		fprintf(tikzfp, "\n};\n");
	}
	fprintf(tikzfp, "\\end{axis}\n\\end{tikzpicture}\n");
	fclose(tikzfp);
}

/**
 * Main function
 */
int main(int argc, char** argv) {
	
	int d = 2;
	int n = 100;
	double sigma = 0.4;
	
	double bmin = 0.0;
	double bmax = -0.22;
	int nb = 30;
	int ns = 5; //Number of different simulations.
	
	minev* mv = (minev*)calloc(1, sizeof(minev));
	mv->ns = ns;
	mv->nb = nb;
	mv->bmin = bmin;
	mv->bmax = bmax;
	
	mv->data = (double*)calloc(nb*(ns+1), sizeof(double));
	mv->bg = (BallGrid*)calloc(ns, sizeof(BallGrid));
	//init_ball_grid(mv->bg+0, d, n, sigma, 3);  //Initializes with different number of grid points.
	init_ball_grid(mv->bg+0, d, n, sigma,  3);
	init_ball_grid(mv->bg+1, d, n, sigma,  6);
	init_ball_grid(mv->bg+2, d, n, sigma, 12);
	init_ball_grid(mv->bg+3, d, n, sigma, 24);
	init_ball_grid(mv->bg+4, d, n, sigma, 48);
	
	run_minev(mv);
	//double* minevls = (double*)calloc(2*nb, sizeof(double));
	//min_eigval_array(bg, bmin, bmax, nb, minevls);
	
	const char fname[] = "clout-test-2d-2";
	//tikz_plot(minevls, nb, 1, "Minimum eigenvalue of $\\matr{A}(\\beta)$.", "$\\beta$~($\\varsigma^{-1}$)", "Minimum eigenvalue", fname);
	tikz_plot_minev(mv, fname);
	
	del_ball_grid(mv->bg+0);
	del_ball_grid(mv->bg+1);
	del_ball_grid(mv->bg+2);
	del_ball_grid(mv->bg+3);
	del_ball_grid(mv->bg+4);
	free(mv->bg);
	free(mv->data);
	free(mv);
	return 0;
}
