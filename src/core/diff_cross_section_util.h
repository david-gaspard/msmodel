/****
 * @date Created on 2021-07-10 at 13:03:41 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the differential cross section utilities.
 ***/
#ifndef _DIFF_CROSS_SECTION_UTIL_H
#define _DIFF_CROSS_SECTION_UTIL_H
#include "diff_cross_section_type.h"
#include "medium_type.h"

void del_diff_cross_section(DiffCrossSection* dcs);
void print_param_diff_cross_section(DiffCrossSection* dcs);
double get_arg_dcs(DiffCrossSection* dcs, int i);
void export_tikz_diff_cross_section(DiffCrossSection* dcs, Medium* med, int thlog, const char* fname);
void save_diff_cross_section(DiffCrossSection* dcs, const char* fname);
void parse_diff_cross_section(DiffCrossSection* dcs, int narg, char** args);

#endif
