/***
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created on 2020-07-25 at 17:59 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Header file of the Bessel and Green functions
 **/
#ifndef _GREEN_BESSEL_H
#define _GREEN_BESSEL_H
#include "dcomplex_type.h"

double uball_volume(int d);
double uball_surface(int d);
dcomplex bessel_k_int(int nu, dcomplex z);
dcomplex bessel_k(int twonu, dcomplex z);
dcomplex hypergeom_0f1_reg_series(double a, dcomplex z);
dcomplex free_green(int d, dcomplex k, double r);
dcomplex free_green_imag(int d, dcomplex k, double r);
dcomplex dk_free_green(int d, dcomplex k, double r);
dcomplex bessel_j_ratio(int twonu, dcomplex z);
dcomplex bessel_k_ratio(int twonu, dcomplex z);

#endif


