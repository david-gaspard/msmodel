/****
 * @date Created on 2021-05-09 at 19:55:53 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the assertion utilities used in tests. Most of these function return 1 on failure, 0 on success.
 ***/
#ifndef _ASSERTION_H
#define _ASSERTION_H
#include "dcomplex_type.h"

void print_test(const char* name);
int assert_vector_equals(int n, double* x, int n_expc, double* x_expc);
int assert_cvector_equals(int n, dcomplex* z, int n_expc, dcomplex* z_expc);
int assert_vector_close(int n, double* x, int n_expc, double* x_expc, double toler);
int assert_cvector_close(int n, dcomplex* z, int n_expc, dcomplex* z_expc, double toler);
int assert_cvector_subset(int n, dcomplex* z, int n_expc, dcomplex* z_expc, double toler);

#endif
