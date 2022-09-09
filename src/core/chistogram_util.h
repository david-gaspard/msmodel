/****
 * @date Created on 2021-09-14 at 20:15:00 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex histogram utilities
 ***/
#ifndef _COMPLEX_HISTOGRAM_UTIL_H
#define _COMPLEX_HISTOGRAM_UTIL_H
#include "medium_type.h"
#include "chistogram_type.h"

void del_chistogram(Chistogram* chist);
int fill_chist(Chistogram* chist, int nz, dcomplex* z);
void print_param_chist(Chistogram* chist);
void export_tikz_chistogram(Chistogram* chist, Medium* med, int id, const char* fname);
void save_chistogram(Chistogram* chist, const char* fname);
void parse_chistogram(Chistogram* chist, int narg, char** args);

#endif
