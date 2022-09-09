/****
 * @date Created on 2020-10-02 at 19:06:31 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Pseudo-random number generator of double precision values in uniform deviate.
 * This generator is in principle thread safe as far as each thread use its own "random" structure.
 * Adapted in C from the book W. H. Press, et al., Numerical Recipes, 3rd ed.
 * Implementation of the highest quality recommended pseudo-random number generator
 * for double precision variables based on 64-bits positive integers. Seed zero is supported.
 * The period of the generator is about 3.138 x 10^57 (much larger than 2^64 or 10^20).
 ***/
#include <math.h>             /* Standard Library for Mathematical Functions */
#include "nru_random_type.h"  /* Import the Internal State of the Numerical Recipes Uniform Generator */

/**
 * Affects the internal state of the random generator, and returns a 64-bit random integer in uniform deviate.
 * See Numerical Recipes book chapter 7, and especially p.342, for explanations of the method.
 */
uint64_t next_state(Urandom* rnd) {
	rnd->u = rnd->u * 2862933555777941757LL + 7046029254386353087LL; //u=Linear Congruential Generator (C3)
	rnd->v ^= rnd->v >> 17; //v=Xorshift Method (A3)
	rnd->v ^= rnd->v << 31;
	rnd->v ^= rnd->v >> 8;
	rnd->w = 4294957665U*(rnd->w & 0xffffffff) + (rnd->w >> 32); //w=Multiply with Carry Method (B1).
	uint64_t x = rnd->u ^ (rnd->u << 21); //x=Another Xorshift Method (A1)
	x ^= x >> 35;
	x ^= x << 4;
	return (x + rnd->v) ^ rnd->w; //Final combination
}

/**
 * Returns a random double-precision floating value in the range 0 to 1 (excluding bounds).
 */
double random_double(Urandom* rnd) {
	//printf("[INFO] RNG Internal state: u = %llu, v = %llu, w = %llu\n", rng->u, rng->v, rng->w); //Just in case of checking.
	return 5.42101086242752217E-20*next_state(rnd); //Just the uint64 random value divided by 2^64.
}

/**
 * Returns a random value according to a normal distribution with standard parameters mu=0, sigma=1.
 * Uses the ratio-of-uniforms rejection-based method (see Numerical Recipes, 3rd ed., p.367).
 */
double random_normal(Urandom* rnd) {
	double u, v, x, y, q;
	do {
		u = random_double(rnd);
		v = 1.7156*(random_double(rnd)-0.5);
		x = u - 0.449871;
		y = fabs(v) + 0.386595;
		q = x*x + y*(0.19600*y-0.25472*x);
	} while (q > 0.27597 && (q > 0.27846 || v*v > -4.*u*u*log(u)));
	return v/u;
}

/**
 * Initialize the pseudo-random generator. Call with any integer seed (except value of v).
 */
void init_random(Urandom* rnd, uint64_t seed) {
	rnd->seed = seed;  //Store the seed.
	rnd->v = 4101842887655102017LL;
	rnd->w = 1;
	rnd->u = seed ^ rnd->v;
	next_state(rnd);
	rnd->v = rnd->u;
	next_state(rnd);
	rnd->w = rnd->v;
	next_state(rnd);
}
