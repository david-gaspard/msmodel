/***
 * @date Created on 2020-08-08 at 11:47 CEST
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Header file providing the utilities related to the random medium.
 **/
#ifndef _MEDIUM_UTIL_H
#define _MEDIUM_UTIL_H
#include "dcomplex_type.h"
#include "medium_type.h"

void init_copy_medium(Medium* srcmed, Medium* dstmed);
void del_medium(Medium* med);
double radius_of_ball(Medium* med);
void fill_medium(Medium* med, uint64_t seed);
char* shape_to_string(Shape shape);
double atom_distance(Medium* med, int i, int j);
double max_distance(Medium* med);
double kimax_value(Medium* med, uint64_t seed);
void center_medium(Medium* med, double* center);
dcomplex invf(Medium* med, dcomplex k);
void print_param_medium(Medium* med);
void print_cross_section(Medium* med, dcomplex k);
dcomplex cross_section(Medium* med, dcomplex k);
dcomplex max_cross_section(Medium* med, dcomplex k);
double mean_free_path(Medium* med, double k);
void detect_roundoff_error(Medium* med, double imk, int verbose);
void export_tikz_medium(Medium* med, const char* fname, uint64_t seed);
void save_medium(Medium* med, const char* fname);
void save_medium_points(Medium* med);
void parse_medium(Medium* med, int narg, char** args);

#endif

