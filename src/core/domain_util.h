/****
 * @date Created on 2021-03-12 at 17:00:40 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex domain utilities.
 ***/
#ifndef _DOMAIN_UTIL_H
#define _DOMAIN_UTIL_H
#include "domain_type.h"

int parse_domain(Domain* dom, int narg, char** args);
void copy_domain(Domain* src, Domain* dst);
void save_domain(Domain* dom, FILE* fp);
double length_x(Domain* dom);
double length_y(Domain* dom);
double domain_area(Domain* dom);
int is_point_in_domain(Domain* dom, double* p);
int is_in_domain(Domain* dom, dcomplex z);
int exclude_from_domain(Domain* dom, int* n, dcomplex* z);
dcomplex halton_point_2_3(Domain* dom, int i);
void set_halton_points(Domain* dom, int n, dcomplex* z);

#endif
