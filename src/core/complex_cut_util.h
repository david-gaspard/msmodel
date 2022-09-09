/****
 * @date Created on 2020-11-30 at 17:06:38 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex cut utilities.
 ***/
#ifndef _COMPLEX_CUT_UTIL_H
#define _COMPLEX_CUT_UTIL_H
#include "complex_cut_type.h"
#include "medium_type.h"

void init_complexcut(ComplexCut* ccut, int ns, dcomplex zmin, dcomplex zmax, int nseed, const char* title);
void del_complexcut(ComplexCut* ccut);
dcomplex get_ccut_arg(ComplexCut* ccut, int l);
dcomplex get_ccut_value(ComplexCut* ccut, int i);
dcomplex get_ccut_laplacian(ComplexCut* ccut, int i);
void print_param_complexcut(ComplexCut* ccut);
void save_complexcut(ComplexCut* ccut, const char* fname);
void export_tikz_ccut(ComplexCut* ccut, Medium* med, int lflag, const char* fname);
void parse_complexcut(ComplexCut* ccut, int narg, char** args);

#endif
