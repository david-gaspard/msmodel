/****
 * @date Created on 2021-08-24 at 13:19:54 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the incomplete gamma function and related functions.
 ***/
#ifndef _INCOMPLETE_GAMMA_H
#define _INCOMPLETE_GAMMA_H

double hypergeom_1f1_series(double a, double b, double z);
double lower_gamma(double a, double z);

#endif
