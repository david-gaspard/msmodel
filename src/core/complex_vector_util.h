/****
 * @date Created on 2021-04-06 at 12:43:52 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex vector utilities.
 ***/
#ifndef _COMPLEX_VECTOR_UTIL_H
#define _COMPLEX_VECTOR_UTIL_H
#include "dcomplex_type.h"

void copy_cvector(int n, dcomplex* src, dcomplex* dst);
double cnorm(int n, dcomplex* x);
dcomplex cnormalize(int n, dcomplex* x);
double cdistance(int n, dcomplex* x, dcomplex* y);
void constant_unit_cvector(int n, dcomplex* x);
dcomplex cdot_product(int n, dcomplex* x, dcomplex* y);
void swap_cvector(dcomplex* z, int i, int j);
void sort_cvector(int n, dcomplex* z, double (*f)(dcomplex));
int uniq_cvector(int* n, dcomplex* z, double toler);
int filter_cvector(int* n, dcomplex* z, int* flag);
void print_cvector(int n, dcomplex* x, const char* name);
int count_complex_data(const char* arg);
void parse_complex_data(int n, dcomplex* data, char* arg);

#endif
