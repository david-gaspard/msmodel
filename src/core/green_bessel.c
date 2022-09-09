/**************************************************************************************************
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created at 2020-07-07 at 17:25 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Computation of the free Green function in arbitrary dimension (typically d<15),
 * and computation of the Bessel function K_nu(z) of small integer orders (nu=0 or around |nu|<6).
 * Uses asymptotic expansion, power series, and the Miller recurrence method to compute the Bessel K function.
 * @ref Inspired by the following references by order of importance:
 * - N. M. Temme, J. Comput. Phys. 19, 326-337 (1975);
 * - I. J. Thompson, A. R. Barnett, Comput. Phys. Commun. 47, 245-257 (1987);
 * - W. H. Press, et al., Numerical Recipes, 3rd ed. (2007);
 * - F. W. J. Olver, et al., NIST Handbook of Mathematical Functions (2010);
 *************************************************************************************************/
//#include <stdio.h>
#include <math.h>       /* Standard Library of Mathematical Functions */
#include "dcomplex_type.h"
/**
 * Define some preprocessor constants. Do not edit unless you know what you are doing.
 */
#define PI        3.1415926535897932385   /* Value of pi constant. */
#define TWOPI     6.2831853071795864769   /* Value of 2*pi. */
#define FOURPI    12.566370614359172954   /* Value of 4*pi. */
#define SQRTPI    1.7724538509055160273   /* Square root of pi. */
#define PSI1     -0.5772156649015328606   /* Digamma function at z=1, i.e., psi(1). Also minus the Euler-Mascheroni gamma constant. */
#define MEPS      1.0e-16                 /* Approximate machine epsilon of double precision. (not used) */
#define MEPS2     1.0e-32                 /* Square of the machine epsilon ((10^-16)^2=10^-32). Used in stopping criteria. */
#define LOGMEPS  -16.0                    /* Base-10 logarithm of the machine epsilon. */
#define LOGMEPS2  256.0                   /* Square of the base-10 logarithm of the machine epsilon (16^2=256). */
#define LNEPS    -36.841361487904730944   /* Neperian logarithm of the machine epsilon (-16*ln(10)). */
#define FOCAL     3.0                     /* Focal length of the parabola of validity of the power series expansion in the complex z-plane. */
#define IMAX      50                      /* Maximum number of terms allowed in the logarithmic series of K_nu(z) at small z (rarely larger than 30). */

/**
 * Computes the volume of the unit d-ball using recursion.
 */
double uball_volume(int d) {
	double vol;
	if (d == 0) {
		return 1.;
	}
	else if (d%2 == 0) {
		vol = PI;
	}
	else {
		vol = 2.;
	}
	int drec = d;
	while (drec > 2) {
		vol  *= TWOPI/drec;
		drec -= 2;
	}
	return vol;
}

/**
 * Computes the surface area of the unit d-ball using recursion.
 */
double uball_surface(int d) {
	double surf;
	if (d == 0) {
		return 0.;
	}
	else if (d%2 == 0) {
		surf = TWOPI;
	}
	else {
		surf = 2.;
	}
	int drec = d;
	while (drec > 2) {
		drec -= 2;
		surf *= TWOPI/drec;
	}
	return surf;
}

/************************************
 * PRIVATE CORE NUMERICAL FUNCTIONS
 ***********************************/
/**
 * Computes the Bessel auxiliary function Btilde_nu1(z)=K_(nu1-0.5)(z)/(sqrt(pi/(2*z))*exp(-z))
 * from the exact representation for odd dimensions. The integer index is nu1=nu+1/2=(d-1)/2.
 * The expansion is finite and exact.
 */
static dcomplex bessel_b_asym_1(int nu1, dcomplex z) {
	dcomplex b = 1.0, t = 1.0;
	dcomplex twoz = 2.0*z;
	int i;
	for (i = 1; i < nu1; i++) {
		t *= (nu1 - i)*(nu1 + i - 1)/(twoz*i);
		b += t;
	}
	//printf("[INFO] B1(nu1=%d, z=%e %+ei) = %.15e %+.15ei\n", nu1, creal(z), cimag(z), creal(b), cimag(b));
	return b;
}

/**
 * Computes the Bessel auxiliary function B_nu(z)=K_nu(z)/(sqrt(pi/(2*z))*exp(-z)) from the asymptotic expansion for integer "nu", i.e., even dimensions (nu=(d-2)/2).
 * The integer order "nu" must be moderate, typically nu<6. This works in double precision for |z|>16.
 */
static dcomplex bessel_b_asym_2(int nu, dcomplex z) {
	dcomplex b = 1.0, t = 1.0;
	dcomplex twoz = 2.0*z;
	int i, imax = (int)floor(cabs(twoz));  //Index at which the term "t" is minimum, the series diverges beyond this.
	for (i = 1; i < imax; i++) {
		t *= (nu - i + 0.5)*(nu + i - 0.5)/(twoz*i);
		b += t;
		if (creal(t)*creal(t) + cimag(t)*cimag(t) < MEPS2) {
			break;
		}
	}
	//printf("[INFO] B2(nu=%d, z=%e %+ei) = %.15e %+.15ei\n", nu, creal(z), cimag(z), creal(b), cimag(b));
	return b;
}

/**
 * Computes the std Bessel function K_nu(z) from the power series at small argument.
 * The region of validity has a parabolic shape with the focus at the origin z=0 and oriented in the direction of the branch cut of K_nu(z).
 */
static dcomplex bessel_k_series(int nu, dcomplex z) {
	int i;
	dcomplex t  = 0.5; //Value for nu=0.
	dcomplex kf = 0.0;
	dcomplex z2 = z*z, z2d4 = z2/4;
	double psinu = PSI1;
	if (nu != 0) {//Prepare for the finite sum in the case nu!=0.
		dcomplex z2rev; //= (-1)^(nu+1)/z^2.
		if (nu%2 == 0) {
			t = 1.0/nu;
			z2rev = -1/z2;
			kf += z2rev/t;
		}
		else {
			t = 1/z;
			z2rev = 1/z2;
		}
		if (nu%4 > 1) {
			kf = -kf;
			t = -t;
		}
		kf += t;
		int imin = nu/2; //Floor integer division.
		for (i = imin; i < nu-1; i++) {//Run half the finite sum.
			t *= z2d4/((i+1)*(i-nu+1));
			kf += t + z2rev/t;
		}
		t *= -z2d4/nu; //Recycle the value of t to feed the logarithmic series (nu!=0).
		//printf("[INFO] Kf(nu=%d, z=%e %+ei) = %.15e %+.15ei\n", nu, creal(z), cimag(z), creal(kf), cimag(kf));
		for (i = 1; i <= nu; i++) {//Computes psi(nu+1).
			psinu += 1.0/i;
		}
		//printf("[INFO] psi(%d) = %.16f\n", nu+1, psinu);
	}
	dcomplex c = PSI1 + psinu - 2*clog(z/2);
	dcomplex kl = c*t;
	//printf("[INFO] t0 = %.15e %+.15ei, kl0 = %.15e %+.15ei\n", creal(t), cimag(t), creal(kl), cimag(kl));
	dcomplex klt;
	for (i = 1; i < IMAX; i++) {//Logarithmic series.
		t *= z2d4/(i*(i + nu));
		c += 1.0/i + 1.0/(i + nu);
		klt = c*t;
		kl += klt;
		if (creal(klt)*creal(klt) + cimag(klt)*cimag(klt) < MEPS2*(creal(kl)*creal(kl) + cimag(kl)*cimag(kl))) {
			//printf("[INFO] Series break at i=%d, with |klt| = %e\n", i, cabs(klt));
			break;
		}
	}
	//if (i == IMAX) {
	//	printf("[WARN] The maximum number of terms (imax=%d) has been reached in the logarihtmic series ! Revise the estimate of the number of terms...\n", IMAX);
	//}
	//printf("[INFO] Kl(nu=%d, z=%e %+ei) = %.15e %+.15ei\n", nu, creal(z), cimag(z), creal(kl), cimag(kl));
	//printf("[INFO] K (nu=%d, z=%e %+ei) = %.15e %+.15ei\n", nu, creal(z), cimag(z), creal(kf + kl), cimag(kf + kl));
	return kf + kl;
}

/**
 * Computes the Bessel auxiliary function B_nu(z)=K_nu(z)/(sqrt(pi/(2*z))*exp(-z)) with the Miller recurrence method [see Temme (1975)].
 * The region of validity is complementary to the power series expansion, and the estimate of the number of iterations only works for |z|<16.
 */
static dcomplex bessel_b_miller(int nu, dcomplex z) {
	dcomplex f = 1., r = 0.;
	int imax = (int)floor((creal(z) - LNEPS)*(creal(z) - LNEPS)/(4*(creal(z) + cabs(z))));
	//printf("[INFO] Miller imax = %d\n", imax);
	double a;
	int i;
	for (i = imax; i > 0; i--) {//Backward recurrence.
		a = i/((i - nu - 0.5)*(i + nu - 0.5));
		f = a*(2*(z + i) - (i + 1)/f);
		r = 1. + r/f;
	}
	//printf("[INFO] BM(nu=%d, z=%e %+ei) = %.15e %+.15ei\n", nu, creal(z), cimag(z), creal(1/r), cimag(1/r));
	return 1./r;
}

/******************************
 * STANDARD BESSEL FUNCTIONS
 *****************************/
/**
 * Computes the Bessel K_nu(z) function for complex "z" and small positive integer "nu" (nu < 6).
 * Choose between the three methods, asymptotic (A), power series (B), continued fraction (C), according to the region of "z"
 * in order to enhance performance and accuracy in double precision arithmetics.
 */
dcomplex bessel_k_int(int nu, dcomplex z) {
	if (creal(z)*creal(z) + cimag(z)*cimag(z) > LOGMEPS2) {//If |z| large enough, then uses the asymptotic expansion.
		return csqrt(0.5*PI/z)*cexp(-z)*bessel_b_asym_2(nu, z);
	}
	else if (creal(z) < FOCAL - cimag(z)*cimag(z)/(4*FOCAL)) {//If Re(z) small enough, then uses the power series expansion.
		return bessel_k_series(nu, z);
	}
	else {//Otherwise uses the Miller recurrence method.
		return csqrt(0.5*PI/z)*cexp(-z)*bessel_b_miller(nu, z);
	}
}

/**
 * Computes the modified Bessel function, K_nu(z), for complex "z" and any integer or half-integer order "nu" passed by "twonu" = 2*nu.
 * Note that the integer value of 2*nu is passed in argument, not "nu" itself.
 */
dcomplex bessel_k(int twonu, dcomplex z) {
	if (twonu < 0) {//Takes the absolute value of "nu" since K_-nu(z) = K_nu(z) for all "nu".
		twonu = -twonu;
	}
	if (twonu == 0 || twonu == 2) {//If nu=0 or nu=1.
		return bessel_k_int(twonu/2, z);
	}
	else if (twonu%2 == 0) {//Integer value of twonu >= 4, i.e., nu >= 2, then use stable forward recurrence.
		dcomplex k0 = bessel_k_int(0, z); //Compute K_0(z).
		dcomplex k1 = bessel_k_int(1, z); //Compute K_1(z).
		dcomplex k2;
		int n, nu = twonu/2;
		for (n = 1; n < nu; n++) {//Stable forward recurrence.
			k2 = (2*n*k1)/z + k0;
			k0 = k1;
			k1 = k2;
		}
		return k2;
	}
	else {//Half-integer value of "nu".
		return csqrt(PI/(2*z))*cexp(-z)*bessel_b_asym_1((twonu + 1)/2, z);
	}
}

/**
 * Computes the regularized hypergeometric function, 0F1(a,z) = Sum(z^n/(n!*gamma(a+n)), for n in[0, inf]), from its power series representation around z=0.
 * The series expansion is accurate in double precision for |z| < 2, roughly speaking. It is also often valid much farther.
 * If "a" is an integer smaller than 1, then use the reflection formula: 0F1(a,z) = z^(1-a) 0F1(2-a,z).
 * Note that the Bessel function is given by J_nu(z) = (z/2)^nu 0F1(nu+1, -(z/2)^2), based on this function 0F1(a,z).
 */
dcomplex hypergeom_0f1_reg_series(double a, dcomplex z) {
	int i;
	dcomplex t = 1., hypf;
	if (a < 1. && fmod(a, 1.) == 0.) {//If integer a<1, use 0F1(a,z) = z^(1-a) 0F1(2-a,z).
		t = cpow(z, 1. - a);
		a = 2. - a;
	}
	t /= tgamma(a); //true gamma function for real arguments.
	hypf = t;
	for (i = 1; i <= IMAX; i++) {
		t *= z/(i*a);
		hypf += t;
		a += 1.;
		if (cabs(t) < MEPS*cabs(hypf)) {
			break;
		}
	}
	return hypf;
}

/************************************
 * GREEN FUNCTIONS AND PROPAGATORS:
 ***********************************/
/**
 * Computes the free Green function for any integer dimension (typically |d|<12) and complex wave number k.
 */
dcomplex free_green(int d, dcomplex k, double r) {
	dcomplex ik = k*I; //Maybe externalize this variable with G(d,ik,r) to avoid recomputation ?
	if (d == 1) {
		return cexp(ik*r)/(2*ik);
	}
	else if (d == 3) {
		return -cexp(ik*r)/(FOURPI*r);
	}
	else if (d == 2) {
		return -bessel_k_int(0, -ik*r)/TWOPI;
	}
	return -cpow(-ik/(TWOPI*r), (d - 2.)/2.)*bessel_k(d - 2, -ik*r)/TWOPI; //Other cases.
}

/**
 * Computes the formal imaginary part of the Green function based on free_green() for complex wave number k.
 * This function is the rigorous analytic continuation of I_k(r) = -Im[G_k(r)] for complex k.
 */
dcomplex free_green_imag(int d, dcomplex k, double r) {
	if (cabs(k)*r < FOCAL) {//For small arguments (FOCAL=3), then use series expansion for to avoid loss of significance.
		dcomplex x = k*r;
		return 0.25*cpow(k/(2*SQRTPI), d-2.)*hypergeom_0f1_reg_series(d/2., -x*x/4);
	}
	int sgn = (creal(k) < 0 && d%2 == 0) ? -1 : 1;   //Sign which is -1 only for Re(k)<0 in even dimensions.
	return sgn*(free_green(d, -k, r) - free_green(d, k, r))/(2*I);
}

/**
 * Derivative with respect to the wave number "k" of the free Green function.
 * This derivative is just given by (k/2pi)*G(d-2, k, r), from the original free Green function G(d,k,r).
 */
dcomplex dk_free_green(int d, dcomplex k, double r) {
	return (k/TWOPI)*free_green(d-2, k, r);
}

/************************************
 * BESSEL FUNCTIONS AND THEIR RATIO:
 ***********************************/
#define REZMAX    50.         /* Maximum value of |Re(z)| beyond which the asymptotic expansion for J_(nu+1)(z)/J_nu(z) will be used, so |Re(z)|>50. */
#define CMAX      200         /* Maximum number of iterations used in continued fractions. */
#define AMAX      200         /* Maximum number of iterations used in the asymptotic expansion of B_nu(z), should be enough for nu<100 and |z|>50. */
/**
 * Computes the Bessel auxiliary function B_nu(z)=K_nu(z)/(sqrt(pi/(2*z))*exp(-z)) from the asymptotic expansion for integer or half-integer "nu".
 * The order "nu" can be large, but not too much, typically nu<100. This expansion is finite and exact for half-integer "nu".
 */
static dcomplex bessel_b_asym(int twonu, dcomplex z) {
	dcomplex t = 1., b = 1., eightz = 8.*z;
	int i, k;
	for (i = 0; i < AMAX; i++) {
		k = 2*i + 1;
		t *= (twonu - k)*(twonu + k)/(eightz*(i+1));
		b += t;
		if (creal(t)*creal(t) + cimag(t)*cimag(t) < MEPS2) {
			break;
		}
	}
	return b;
} 

/**
 * Computes the ratio of Bessel functions, J_(nu+1)(z)/J_nu(z), for integer or half-integer "nu" and complex "z" using the asymptotic expansion.
 * The variable "z" must be |z|>16 for this method to work in double precision.
 */
static dcomplex bessel_j_ratio_asym(int twonu, dcomplex z) {
	int sgn = 1;
	if (creal(z) < 0) {//The function is odd, and the left half plane involves a different formula.
		z = -z;
		sgn = -1;
	}
	dcomplex iph = I*(z - PI*(twonu - 1.)/4.);
	dcomplex epi = cexp(iph), emi = 1./epi;
	return sgn*I*(emi*bessel_b_asym(twonu+2, I*z) + epi*bessel_b_asym(twonu+2, -I*z))/(emi*bessel_b_asym(twonu, I*z) - epi*bessel_b_asym(twonu, -I*z));
}

/**
 * Computes the ratio of Bessel functions, J_(nu+1)(z)/J_nu(z), for integer or half-integer "nu" and complex "z" using continued fraction.
 * The continued fraction is evaluated with the Lentz method (see Numerical Recipes). The argument |z| must be small enough compared to CMAX.
 * In fact, this algorithm extends to any complex "nu", but we do not need such generalization here.
 */
static dcomplex bessel_j_ratio_cfrac(int twonu, dcomplex z) {
	dcomplex fac, f = MEPS2, c = f, d = 0.;
	dcomplex a = -z*z;
	dcomplex b = twonu + 2.;
	int i;
	for (i = 0; i < CMAX; i++) {//Continued fraction evaluation using the Lentz method.
		b += 2.;
		d = b + a*d;
		if (d == 0.)
			d = MEPS2;
		d = 1./d;
		c = b + a/c;
		if (c == 0.)
			c = MEPS2;
		fac = c*d;
		f *= fac;
		if (creal(fac - 1.)*creal(fac - 1.) + cimag(fac)*cimag(fac) < MEPS2)
			break;
	}
	//if (i >= CMAX/2)
	//	printf("[WARN] Exceeding iteration i=%d, for J(nu=%g, z=%g%+gi).\n", i, twonu/2., creal(z), cimag(z));
	return z/(twonu + 2. + f);
}

/**
 * Computes the ratio of Bessel functions, J_(nu+1)(z)/J_nu(z), using either the continued fraction for small "z" or the asymptotic expansion.
 */
dcomplex bessel_j_ratio(int twonu, dcomplex z) {
	if (fabs(creal(z)) > REZMAX) {//If |Re(z)| too large, then uses the asymptotic expansion.
		return bessel_j_ratio_asym(twonu, z);
	}
	else {//Otherwise, use the continued fraction.
		return bessel_j_ratio_cfrac(twonu, z);
	}
}

/**
 * Computes the ratio of Bessel functions, K_(nu+1)(z)/K_nu(z), for integer or half-integer "nu" and complex "z" using continued fraction.
 * The continued fraction is evaluated with the Lentz method. The variable "z" must be far enough from the branch cut lying on the negative real axis.
 */
static dcomplex bessel_k_ratio_cfrac(int twonu, dcomplex z) {
	dcomplex fac, f = MEPS2, c = f, d = 0.;
	dcomplex a = (twonu - 1.)*(twonu + 1.);
	dcomplex b = 4.*z;
	int i;
	for (i = 0; i < CMAX; i++) {//Continued fraction evaluation using the Lentz method.
		a -= 8.*i;
		b += 4.;
		d = b + a*d;
		if (d == 0.)
			d = MEPS2;
		d = 1./d;
		c = b + a/c;
		if (c == 0.)
			c = MEPS2;
		fac = c*d;
		f *= fac;
		if (creal(fac - 1.)*creal(fac - 1.) + cimag(fac)*cimag(fac) < MEPS2)
			break;
	}
	//if (i >= CMAX/2)
	//	printf("[WARN] Exceeding iteration i=%d, for K(nu=%g, z=%g%+gi).\n", i, twonu/2., creal(z), cimag(z));
	return 1. + (twonu + 1.)/(2*z) + f/(2*z);
}

/**
 * Computes the ratio of modified Bessel functions, K_(nu+1)(z)/K_nu(z), using either the continued fraction
 * or the series expansion for small "z" near the branch cut.
 * Uses the fact that K_-nu(z) = K_nu(z) to make "nu" positive.
 */
dcomplex bessel_k_ratio(int twonu, dcomplex z) {
	int refl = (twonu < -1); //Reflection on "nu" must take place.
	if (refl) {//Reflection formula for "twonu".
		twonu = -twonu - 2;
	}
	dcomplex kratio;
	if (twonu%2 == 0 && creal(z) < FOCAL - cimag(z)*cimag(z)/(4*FOCAL) && creal(z) > -REZMAX) {//If z is close to the branch cut and integer nu, then use the power series !
		kratio = bessel_k_series((twonu+2)/2, z)/bessel_k_series(twonu/2, z);
	}
	else {//Otherwise, use the continued fraction.
		kratio = bessel_k_ratio_cfrac(twonu, z);
	}
	return refl ? 1./kratio : kratio;
}
