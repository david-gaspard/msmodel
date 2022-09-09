/****
 * @date Created on 2021-07-31 at 19:55:47 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the utilities to manipulate the histogrma of Im(mu).
 ***/
#ifndef _IMAG_MU_HISTOGRAM_UTIL_H
#define _IMAG_MU_HISTOGRAM_UTIL_H
#include "imag_mu_histogram_type.h"
#include "medium_type.h"
#include "dcomplex_type.h"

void del_ihist(ImagMuHistogram* ihist);
void print_param_ihist(ImagMuHistogram* ihist);
void setup_domain_ihist(ImagMuHistogram* ihist, Medium* med);
int append_list_ihist(ImagMuHistogram* ihist, int nz, dcomplex* z);
void normalize_ihist(ImagMuHistogram* ihist);
void export_tikz_ihist(ImagMuHistogram* ihist, Medium* med, const char* fname);
void save_ihist(ImagMuHistogram* ihist, const char* fname);
void parse_ihist(ImagMuHistogram* ihist, int narg, char** args);

#endif
