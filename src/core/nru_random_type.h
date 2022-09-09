/****
 * @date Created on 2020-10-02 at 19:25:57 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the internal state of the Numerical Recipes highest quality generator.
 ***/
#ifndef _NRU_RANDOM_TYPE_H
#define _NRU_RANDOM_TYPE_H
#include <stdint.h> /* Standard Library for Fixed-Size Integers */

typedef struct NRURandomState_s {//Internal state of the Numerical Recipes Uniform Generator.
	uint64_t u;    //State variable of the Linear Congruential Generator (LCG).
	uint64_t v;    //State variable of the Xorshift method (XOR).
	uint64_t w;    //State variable of the Multiply with Carry method (MWC).
	uint64_t seed; //Initial seed of the generator. This value is not altered state after state.
} Urandom;

#endif
