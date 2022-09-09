/****
 * @date Created on 2021-02-23 at 16:54:15 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the complex cut utilities with logarithmic grid points.
 ***/
#include <stdlib.h>             /* Standard Library for Memory Allocation */
#include <stdio.h>              /* Standard Library for Input and Output */
#include <string.h>             /* Standard Library for String Manipulation */
#include <math.h>               /* Standard Library for Mathematical Functions */
#include "common_util.h"        /* Import General Purpose Functions such as the Progress Bar */
#include "log_cut_type.h"       /* Import the Complex Cut Type */
#include "incomplete_gamma.h"   /* Import the Incomplete Gamma and Related functions. */
#include "complex_data_util.h"  /* Import the Complex Data Parsing functions */
#include "medium_util.h"        /* Import the Medium Utilities, used in the plotting functions */
#include "classical_diffusion_util.h"  /* Import the Classical Diffusion Solver Functions, used in plot */
#include "tikz_util.h"          /* Import TikZ Plots Functions */
/**
 * Define some constant macros:
 */
#define TWOPI       6.2831853071795864769   //Constant 2*pi.
#define SQRTTWOPI   2.5066282746310005024   //Square root of 2*pi.
/**
 * LAPACK's routine to compute the solution to a real system of linear equations.
 */
void dgesv_(int* n, int* nrhs, double* a, int* lda, int* ipiv, double* b, int* ldb, int* info);
/**
 * Returns the position in the k-plane corresponding to the given index "i" in the range {0, ..., ns-1}.
 */
dcomplex get_arg_lcut(LogCut* lcut, int i) {
	return lcut->k[i];
}
/**
 * Returns the real-valued curvilinear coordinate of the given sample point k[i].
 * Since the log cut is oriented in the imaginary direction, this is just the absolute value of Im(k[i]).
 */
double get_coord_lcut(LogCut* lcut, int i) {
	return fabs(cimag(lcut->k[i]));
}
/**
 * Returns the second derivative of the fast cut along the cut direction at the given index "i" in the range {0, ..., ns-1}.
 * This is only an approximation of the true Laplacian restricted to functions which are nearly constant/linear in the transverse direction.
 * Resets the index "i" in the interior of the grid, i.e., at least i=1 and at most i=(ns-2).
 * Calls LAPACK's dgesv_() subroutine to locally fit the unequally spaced data with a parabola, then extracts the second derivative.
 */
double get_laplacian_lcut(LogCut* lcut, int i) {
	if (i <= 0) {
		i = 1;
	}
	else if (i >= lcut->ntot-1) {
		i = lcut->ntot-2;
	}
	double x1 = get_coord_lcut(lcut, i-1); //Calculates the curvilinear coordinates of the three needed points.
	double x2 = get_coord_lcut(lcut, i);
	double x3 = get_coord_lcut(lcut, i+1);
	double a[] = {1., 1., 1., x1, x2, x3, x1*x1/2, x2*x2/2, x3*x3/2};  //Vandermonde matrix.
	double b[] = {creal(lcut->data[i-1]), creal(lcut->data[i]), creal(lcut->data[i+1])};
	int n = 3, nrhs = 1, ipiv[n], info;
	dgesv_(&n, &nrhs, a, &n, ipiv, b, &n, &info);  //Solves the system.
	return b[2];  //Returns the second derivative.
}

/**
 * Initialization function of an empty logarithmic cut with the given parameters.
 * The content of the fast cut must be deleted with the corresponding del_*() function.
 */
void init_logcut(LogCut* lcut, int ns, double rek, double imkmin, double imkmax, int nseed, const char* title) {
	if (ns > 10000000) {//Beyond 10 MegaPixels is a bit too large...
		printf("[ERROR] The number of samples is too large (ns=%d), aborting...\n", ns);
		exit(EXIT_FAILURE);
	}
	if (ns < 2) {
		printf("[ERROR] The number of samples is too small (ns=%d), aborting...\n", ns);
		exit(EXIT_FAILURE);
	}
	if (imkmin >= 0 || imkmax >= 0) {
		printf("[ERROR] The end points must be negative (imkmin=%g, imkmax=%g), aborting...\n", imkmin, imkmax);
		exit(EXIT_FAILURE);
	}
	lcut->rek = rek;
	lcut->imkmin = imkmin;
	lcut->imkmax = imkmax;
	lcut->ntot = ns;
	lcut->nseed = nseed;
	strncpy(lcut->title, title, 79); //Maximum size of the title.
	lcut->realtime = 0.0; //Computation time is zero, since empty.
	lcut->k = (dcomplex*)calloc(lcut->ntot, sizeof(dcomplex));
	int i;
	double lfac = log(imkmax/imkmin);
	for (i = 0; i < lcut->ntot; i++) {//Prepares the logarithmic mesh.
		lcut->k[i] = rek + I*imkmin*exp((lfac*i)/(lcut->ntot-1));
	}
	lcut->data = (dcomplex*)calloc(lcut->ntot, sizeof(dcomplex));
	lcut->get_arg = (dcomplex(*)(void*, int))get_arg_lcut;
}

/**
 * Deletes the content (allocated pointers) within the given fast cut, not the pointer itself.
 * Do not forget to invoke after use of "init_*" functions.
 */
void del_logcut(LogCut* lcut) {
	free(lcut->k);
	free(lcut->data);
}

/**
 * Prints a short summary of the parameters used in the logarithmic cut for the standard output.
 */
void print_param_logcut(LogCut* lcut) {
	printf("[INFO] Log cut entitled '%s' with %d samples at rek=%g, imkmin=%g, imkmax=%g, for %d seeds.\n",
		lcut->title, lcut->ntot, lcut->rek, lcut->imkmin, lcut->imkmax, lcut->nseed);
}

/**
 * Saves the log cut data, in addition to other metadata, to the file titled "fname".
 */
void save_logcut(LogCut* lcut, const char* fname) {
	FILE* fp = fopen(fname, "a"); //Append data to the given file.
	fprintf(fp, "[logkcut]\ntitle=%s\n", lcut->title);  //Title of the log cut.
	fprintf(fp, "realtime=%g\n", lcut->realtime); //Total physical duration of the computation of the log cut (in seconds).
	fprintf(fp, "nseed=%d\n", lcut->nseed);
	fprintf(fp, "sample=%d\n", lcut->ntot);
	fprintf(fp, "rek=%.12g\n", lcut->rek);
	fprintf(fp, "imkmin=%.12g\n", lcut->imkmin);
	fprintf(fp, "imkmax=%.12g\n", lcut->imkmax);
	save_cdata((ComplexData*)lcut, fp);  //Call to super-class function to store the data.
	fprintf(fp, "\n");
	fclose(fp);
}

/**
 * Returns the upper bound on the potential function.
 */
double potential_upper_bound(Medium* med, dcomplex k) {
	dcomplex f = invf(med, k);
	double invf2 = creal(f)*creal(f) + cimag(f)*cimag(f); //Square modulus of the inverse scattering amplitude.
	double d = (double)(med->d);  //Dimension of the medium.
	double c = (d + 3.)/2;   //Corner parameter for a d-ball.
	double diam = 2*radius_of_ball(med);  //Diameter of the d-ball.
	double x = -2*cimag(k)*diam; //Argument of the hypergeometric function that follows.
	double hypf = pow(x, -c)*exp(x)*lower_gamma(c+1, x); //Computes 1F1(1, c+1, -2Im(k)L) for Im(k)<0 from the lower incomplete gamma function.
	double avgg2 = (pow(cabs(k), d-3.)/(4.*pow(TWOPI*diam, d-1.)))*(tgamma(c+d)/(tgamma(c+1)*tgamma(d)))*hypf; //Approximate average of |G(k,s)|^2 in a d-ball.
	return 0.5*log(invf2 + (med->n-1)*avgg2);
}

/**
 * Appends a plot in log-log scale of the log cut to the tikz picture file "fname.tikz".
 * Shows both the function and its second derivative.
 */
void export_tikz_lcut(LogCut* lcut, Medium* med, const char* fname) {
	char tikzfile[strlen(fname)+6];
	sprintf(tikzfile, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfile, "a");
	//Initializes the data to be plotted:
	Curve* plot1 = (Curve*)calloc(1, sizeof(Curve)); //Plot1 contains the points to show in the upper panel (resonance density).
	init_curve(plot1, lcut->ntot, "blue,  thick, mark=*, mark size=0.5, line join=bevel");  //Numerical resonance density.
	Curve* plot2 = (Curve*)calloc(2, sizeof(Curve)); //Plot2 contains the points to show in the lower panel (resonance potential).
	init_curve(plot2+0, lcut->ntot, "blue,  thick, mark=*, mark size=0.5");  //Numerical resonance potential function.
	init_curve(plot2+1, lcut->ntot, "black, thick, densely dashed");  //Upper bound of the potential function.
	int i;
	double xi;
	for (i = 0; i < lcut->ntot; i++) {//Loop on the points of the log cut to be plotted.
		xi = get_coord_lcut(lcut, i);  //Curvilinear coordinate xi=-Im(k) (in units of "k").
		//Upper panel of the figure (resonance density):
		plot1[0].x[i] = xi;
		plot1[0].y[i] = get_laplacian_lcut(lcut, i); //Calculate the second derivative (resonance density).
		//Lower panel of the figure (potential function):
		plot2[0].x[i] = xi;
		plot2[0].y[i] = creal(lcut->data[i]); //Extract the real part (potential function).
		plot2[1].x[i] = xi;
		plot2[1].y[i] = potential_upper_bound(med, get_arg_lcut(lcut, i)); //Upper bound on the potential function.
	}
	//Prepares the plot options:
	double xmin = get_coord_lcut(lcut, 0);
	double xmax = get_coord_lcut(lcut, lcut->ntot-1);
	char axisopt1[400]; //Axis options of the upper panel.
	sprintf(axisopt1, "name={reson-density},\n\txmode=log, ymode=log,\n\txmin=%g, xmax=%g,\n\ttitle={%dD %s for $N=%d$ and $k_{\\rm r}=%g\\,\\varsigma^{-1}$\\\\ $\\ell=%.4g\\,\\varsigma$ with %d seeds (%.1f s)},\n\tylabel={$2\\pi\\rho^{(2)}(k)=\\nabla^2\\Phi(k)$},\n\txticklabel=\\empty,",
		xmin, xmax, med->d, shape_to_string(med->shape), med->n, lcut->rek, mean_free_path(med, lcut->rek), lcut->nseed, lcut->realtime);
	char axisopt2[300];  //Axis options of the lower panel.
	sprintf(axisopt2, "name={reson-potential},\n\tat={(reson-density.below south west)},\n\tanchor={above north west},\n\txmode=log,\n\txmin=%g, xmax=%g,\n\tylabel={$\\Phi(k)$},\n\txlabel={$-\\Im k\\varsigma$},", xmin, xmax);
	
	char precmd[900];
	int len = 0; //Length of "precmd".
	if (med->shape == ball || (med->shape == cube && med->d == 1)) {//If the medium is ball shaped, then shows the classical estimate of the peak position.
		double radius = radius_of_ball(med);  //Radius of the d-ball.
		double nsigma = cabs(cross_section(med, lcut->rek));
		len += diffusion_result_plot_v2(med->d, radius, nsigma, precmd);  //Solves the diffusion, and writes the plotting commands to "precmd" (v1=Original version from the draft notes, v2=Version of the final report of the thesis).
	}
	//Computes the range of validity of the potential function method due to round-off errors:
	double abskimax = fabs(kimax_value(med, 1)); //Compute kimax in absolute value with seed=1.
	len += sprintf(precmd+len, "\\draw[black!50] (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymin})}) -- (axis cs:%g,{exp(\\pgfkeysvalueof{/pgfplots/ymax})}) node[pos=0,above right,rotate=90]{$k_{\\rm imax}$}; %%Range of validity of the potential function method.\n", abskimax, abskimax);
	//Adds the approximate power-law behavior given by the upper bound (should only work for ball shaped medium, but approximates everything):
	len += sprintf(precmd+len, "\\addplot[black, thick, densely dashed, domain=\\pgfkeysvalueof{/pgfplots/xmin}:\\pgfkeysvalueof{/pgfplots/xmax}, samples=3] {%g/x^2}; \\addlegendentry{$(d+3)/4k_{\\rm i}^2$}\n", (med->d+3)/4.);
	//else if (med->shape == gaussian) {//Approximate behavior given by the upper bound for a gaussian medium.
	//	double l = pow(med->n, 1./med->d)/SQRTTWOPI; //Computes the standard deviation of the gaussian density around the medium center (r=0).
	//	len += sprintf(precmd+len, "\\addplot[red, thick, domain=\\pgfkeysvalueof{/pgfplots/xmin}:\\pgfkeysvalueof{/pgfplots/xmax}, samples=3] {%g}; \\addlegendentry{$4L^2$}\n", 4*l*l);
	//}
	fprintf(tikzfp, "\\par\\begin{tikzpicture}[%%\n\t/pgfplots/width=250pt,\n\t/pgfplots/height=150pt,\n\t/pgfplots/enlarge x limits=false,\n\t/pgfplots/title style={align=center},\n\t/pgfplots/yticklabel style={rotate=90},\n\t/pgfplots/ticklabel style={%%\n\t\tscaled ticks=false,\n\t\t/pgf/number format/fixed,\n\t\t/pgf/number format/precision=3,\n\t\t/pgf/number format/1000 sep={\\,},\n\t},\n]%%\n");
	tikz_plot(tikzfp, plot1, 1, axisopt1, precmd);  //Plotting the upper panel.
	tikz_plot(tikzfp, plot2, 2, axisopt2, "");      //Plotting the lower panel.
	fprintf(tikzfp, "\\end{tikzpicture}\n");
	
	del_curve(plot1);
	del_curve(plot2+0);
	del_curve(plot2+1);
	free(plot1);
	free(plot2);
	fclose(tikzfp);
}

/**
 * Parses the log cut "lcut" with the arguments of "args".
 * The "lcut" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_logcut(LogCut* lcut, int narg, char** args) {
	int nargmin = 5; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(1);
	}
	int ns, nseed;
	double rek, imkmin, imkmax;
	char *title, *value1;
	value1 = get_value(narg, args, "sample");
	if (sscanf(value1, "%d", &ns) != 1) {
		printf("[ERROR] Invalid number of samples '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "rek");
	if (sscanf(value1, "%lg", &rek) != 1) {
		printf("[ERROR] Invalid real part of k '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "imkmin");
	if (sscanf(value1, "%lg", &imkmin) != 1) {
		printf("[ERROR] Invalid endpoint '%s', aborting...\n", value1);
		exit(1);
	}
	value1 = get_value(narg, args, "imkmax");
	if (sscanf(value1, "%lg", &imkmax) != 1) {
		printf("[ERROR] Invalid endpoint '%s', aborting...\n", value1);
		exit(1);
	}
	
	value1 = get_value(narg, args, "nseed");
	if (sscanf(value1, "%d", &nseed) != 1) {
		printf("[ERROR] Invalid number of seeds '%s', aborting...\n", value1);
		exit(1);
	}
	title = get_value(narg, args, "title");
	init_logcut(lcut, ns, rek, imkmin, imkmax, nseed, title); //Call initialization function.
	sscanf(get_value(narg, args, "realtime"), "%lg", &(lcut->realtime)); //Set realtime after initialization.
	value1 = get_value(narg, args, "data");  //Now parsing the list of data.
	parse_cdata((ComplexData*)lcut, value1);  //Call to super-class function to parse the data string.
}
