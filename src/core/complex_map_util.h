/****
 * @date Created on 2020-08-12 at 11:16:37 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex map management utilities.
 ***/
#ifndef _COMPLEX_MAP_UTIL_H
#define _COMPLEX_MAP_UTIL_H
#include "complex_map_type.h"
#include "medium_type.h"

void init_complexmap(ComplexMap* cmap, int nx, int ny, double xmin, double xmax, double ymin, double ymax, int nseed, const char* title);
void del_complexmap(ComplexMap* cmap);
void find_min_max_cmap(ComplexMap* cmap, int verbose);
dcomplex get_arg_cmap(ComplexMap* cmap, int l);
dcomplex get_value_cmap(ComplexMap* cmap, int ix, int iy);
void set_value_cmap(ComplexMap* cmap, int ix, int iy, dcomplex val);
dcomplex interpolate_bilinear(ComplexMap* cmap, dcomplex z);
void init_laplacemap(ComplexMap* cmap, ComplexMap* laplmap);
void print_param_complexmap(ComplexMap* cmap, int verbose);
void print_complexmap(ComplexMap* cmap, int s); //Deprectated, for test purpose only.
void save_complexmap(ComplexMap* cmap, const char* fname);
void export_tikz_cmap(ComplexMap* cmap, Medium* med, int lmax, int id, const char* fname);
void export_tikz_cut(ComplexMap* cmap, double x, const char* fname);
void parse_complexmap(ComplexMap* cmap, int narg, char** args);

#endif
