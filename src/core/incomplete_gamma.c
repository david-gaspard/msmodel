/****
 * @date Created on 2021-08-24 at 12:50:10 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the incomplete gamma function and related functions for real arguments.
 ***/
#include <stdio.h>
#include <math.h>
/**
 * Define some preprocessor constants. Do not edit unless you know what you are doing.
 */
#define PI        3.1415926535897932385   /* Value of pi constant. */
#define MEPS      1.e-16                  /* Approximate machine epsilon of double precision. */
#define RASYM     18.                     /* Heuristic radius of validity of the asymptotic expansion in double precision, probably related to LNEPS. */
#define IMAX      100                     /* Maximum number of terms allowed in the power series expansions. */

/**
 * Computes the confluent hypergeometric function 1F1(a,b,z) = Sum((a)_n z^n/(b)_n n!), for n in[0, inf]), from its power series representation around z=0.
 * The series expansion is accurate in double precision for |z| < 2, roughly speaking. It is also often valid much farther in "z", especially when a>0, b>0, z>0.
 */
double hypergeom_1f1_series(double a, double b, double z) {//Internal function.
	double t = 1., f = t;
	int i;
	for (i = 1; i <= IMAX; i++) {
		t *= a*z/(b*i);
		f += t;
		a += 1.;
		b += 1.;
		if (fabs(t) < MEPS*fabs(f)) {//Stopping criterion.
			break;
		}
	}
	return f;
}

/**
 * Computes the lower incomplete gamma function using its power series expansion at z=0.
 * This function is typically valid in the interval 0 < z < RASYM, but is only limited by the maximum number of terms IMAX.
 */
double lower_gamma_series(double a, double z) {//Internal function.
	double al = a, t = 1., f = t;
	int i;
	for (i = 0; i < IMAX; i++) {
		al += 1.;
		t *= z/al;
		f += t;
		if (fabs(t) < MEPS*fabs(f)) {//Stopping criterion.
			break;
		}
	}
	return pow(z, a)*exp(-z)*f/a;
}

/**
 * Computes the lower incomplete gamma function using its power series expansion at z=0.
 * This function is only valid for large "z", typically |z| > RASYM. It would also be valid for complex "z".
 */
double lower_gamma_asym(double a, double z) {//Internal function.
	int i, imax = (int)round(fabs(z) + a); //Estimates the index of the smallest term.
	double al = a, t = 1./z, f = t;
	for (i = 0; i < imax; i++) {
		al -= 1.;
		t *= al/z;
		f += t;
		if (fabs(t) < MEPS*fabs(f)) {//Stopping criterion.
			break;
		}
	}
	return tgamma(a) - pow(z, a)*exp(-z)*f;
}

/**
 * Computes the lower incomplete gamma function for positive real arguments "z".
 */
double lower_gamma(double a, double z) {
	if (z < 0.) return NAN; //Complex unimplemented region.
	if (z > RASYM) {//If z is large enough, then use the asymptotic expansion.
		return lower_gamma_asym(a, z);
	}
	return lower_gamma_series(a, z); //Otherwise, use the power series expansion.
}

