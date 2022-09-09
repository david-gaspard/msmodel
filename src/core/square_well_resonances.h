/****
 * @date Created on 2021-03-12 at 10:58:44 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the square well resonances root-finding utilities.
 ***/
#ifndef _SQUARE_WELL_RESONANCES_H
#define _SQUARE_WELL_RESONANCES_H
#include "medium_type.h"
#include "domain_type.h"

void find_resonances(Medium* med, Domain* dom, int lmax, int* nr, dcomplex** reson);
char* resonance_plot_string(Medium* med, Domain* dom, int lmax);

#endif
