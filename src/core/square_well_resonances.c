/****
 * @date Created on 2021-03-11 at 17:41:21 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to find all the complex resonances of the square well approximation of the ball-shaped random Lorentz gas.
 * This approximation model is mostly physically relevant at low energy, i.e., |k|*sp < pi, where "sp" is the mean inter-atomic spacing (unit length).
 * Use Pade approximation to first estimate the roots of the resonance equations, then polish them using the secant method.
 * The resonances are found up to a given maximum angular momentum partial wave "lmax".
 ***/
#include <stdlib.h>         /* Standard Library for Memory Allocation */
#include <stdio.h>          /* Standard Library for Input and Output */
#include <math.h>           /* Standard Library of Mathematical Functions */
#include "green_bessel.h"   /* Import the Library to Compute Ratios of Bessel functions */
#include "domain_util.h"    /* Import the Library to Manipulate Complex Domains */
#include "medium_util.h"    /* Import the Medium Utilities */
#include "root_finder.h"    /* Import the Root-Finding Utilities */

/**
 * Evaluates the resonance function to be root-finded at a given partial wave "l" and at a given complex wave number "k".
 * It can read, i*k*K_(nu+1)(-i*k*R)/K_nu(-i*k*R) + kappa*J_(nu+1)(kappa*R)/J_nu(kappa*R), where nu = l + (d-2)/2,
 * kappa = sqrt(k^2 - F(k)) is the effective wave number in the medium, and R is the radius of the medium.
 */
dcomplex resonance_function(Medium* med, int l, dcomplex k) {
	int twonu = 2*l + med->d - 2;  //Two times nu = l + (d-2)/2.
	dcomplex ik = I*k;
	dcomplex kappa = csqrt(k*k - 1./invf(med, k));  //Effective complex wave number in the medium, since the atomic density is unity.
	double radius = radius_of_ball(med);  //Radius of the medium, assumed to be ball shaped.
	return ik*bessel_k_ratio(twonu, -ik*radius) + kappa*bessel_j_ratio(twonu, kappa*radius); //Returns the result.
}

/**
 * Main function finding the resonances of the effective square well model in a given complex domain "dom" and up to
 * a certain partial wave "lmax" (included). This approximation is mainly physically relevant at low energy, |k|sp < pi.
 */
void find_resonances(Medium* med, Domain* dom, int lmax, int* nr, dcomplex** reson) {
	int maxit = 30;        //Maximum number of root-finding iterations.
	double toler = 1e-12;  //Tolerance on the relative errors of the roots.
	int verb = 0;  //Verbosity level of the root finder (0=quiet, 1=verbose, 2=debug).
	int l, nz;
	dcomplex f(dcomplex k) {
		return resonance_function(med, l, k);
	}
	for (l = 0; l <= lmax; l++) {//Loop on the partial waves "l" to be solved for the resonances.
		nz = nr[l];  //Store the estimated number of resonances since changed by the solve*() routine.
		solve_from_guess(f, dom, &nr[l], reson[l], maxit, toler, verb);
		printf("[INFO] Partial wave l=%d, found %d/%d resonances...\n", l, nr[l], nz);
		fflush(stdout);
	}
}

/**
 * Computes the resonances of the effective square well model up to partial wave "lmax">0 (included), and returns the corresponding TikZ/PGFPlots plot string.
 * The output string "str" should be freed externally after used.
 */
char* resonance_plot_string(Medium* med, Domain* dom, int lmax) {
	int nz = (int)ceil(length_x(dom)*radius_of_ball(med) + 1);  //Roughly estimated number of resonances for each partial wave with safety margin.
	if (med->d%2 == 0) {//In even dimensions, there is a branch cut in the Bessel K_nu(z), so more samples are needed to fit the branch cut with pairs poles/zeros. 
		nz = 3*nz;
	}
	int* nr = (int*)calloc(lmax+1, sizeof(int));  //Array of numbers of valid resonances [nr_0, nr_1, ..., nr_lmax].
	dcomplex** reson = (dcomplex**)calloc(lmax+1, sizeof(dcomplex*));  //List of resonances stored as [z1_0, z2_0, ..., znz_0], [z1_1, ...], ..., [..., znz_lmax].
	int l;
	for (l = 0; l <= lmax; l++) {//Allocate space for the resonances of each partial wave.
		nr[l] = nz;  //Initialize the estimated number of resonances.
		reson[l] = (dcomplex*)calloc(nr[l], sizeof(dcomplex));
	}
	find_resonances(med, dom, lmax, nr, reson);  //Call the resonance root-finder.
	int i, maxlen = (lmax+1)*(200 + 40*nz);  //Estimated output string length with safety margins.
	char* str = (char*)calloc(maxlen, sizeof(char)); //Allocate output string.
	int len = 0;  //Current length of the string.
	len += sprintf(str+len, "\\addplot[mark=o, only marks, mark size=1, nodes near coords, point meta=explicit symbolic] coordinates {%%square well approx resonances up to lmax=%d.\n", lmax);
	for (l = 0; l <= lmax; l++) {//Print the TikZ code to the string for each partial wave.
		for (i = 0; i < nr[l]; i++) {//Loop on confirmed resonances.
			len += sprintf(str+len, "\t(%.12g, %.12g) [%d]\n", creal(reson[l][i]), cimag(reson[l][i]), l);
		}
		free(reson[l]);
	}
	len += sprintf(str+len, "};\n");
	free(reson);
	free(nr);
	return str;
}

/*** OLD PARTS ***
static const int tikzstylelen = 12;
static const char* tikzstyle[] = {
	"mark=o, mark size=1",
	"mark=square, mark size=0.9",
	"mark=triangle, mark size=1.2",
	"mark=diamond, mark size=1.2",
	"mark=pentagon, mark size=1",
	"mark=+, mark size=1.2",
	"mark=star, mark size=1.2",
	"mark=x, mark size=1.3",
	"mark=Mercedes star, mark size=1.2",
	"mark=oplus, mark size=1",
	"mark=10-pointed star, mark size=1.2",
	"mark=otimes, mark size=1",
};
for (l = 0; l <= lmax; l++) {//Print the TikZ code to the string for each partial wave.
	len += sprintf(str+len, "\\addplot[%s] coordinates {%%l=%d, square well approx resonances.\n\t", tikzstyle[l%tikzstylelen], l);
	for (i = 0; i < nr[l]; i++) {//Loop on confirmed resonances.
		len += sprintf(str+len, "(%.12g, %.12g) ", creal(reson[l][i]), cimag(reson[l][i]));
	}
	len += sprintf(str+len, "\n}; \\addlegendentry{$\\ell=%d$}\n", l);
	free(reson[l]);
}
*/
