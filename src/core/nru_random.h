/****
 * @date Created on 2020-10-02 at 19:29:24 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 ***/
#ifndef _NRU_RANDOM_H
#define _NRU_RANDOM_H
#include "nru_random_type.h"  /* Import the Internal State of the Numerical Recipes Uniform Generator */

double random_double(Urandom* rnd);
double random_normal(Urandom* rnd);
void init_random(Urandom* rnd, uint64_t seed);

#endif
