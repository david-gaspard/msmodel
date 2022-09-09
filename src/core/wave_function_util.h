/****
 * @date Created on 2021-11-28 at 12:44:11 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the utilities to manipulate and plot the wave function structure.
 ***/
#ifndef _WAVE_FUNCTION_UTIL_H
#define _WAVE_FUNCTION_UTIL_H
#include "wave_function_type.h"   /* Import the Wave Function Type */
#include "medium_type.h"          /* Import the Medium Type */

void get_pos_wfun(WaveFunction* wfun, int l, double* pos);
void del_wfun(WaveFunction* wfun);
void divide_wfun(WaveFunction* wfun, double div);
void print_param_wfun(WaveFunction* wfun);
void export_tikz_wfun(WaveFunction* wfun, Medium* med, int id, const char* fname);
void save_wfun(WaveFunction* wfun, const char* fname);
void parse_wfun(WaveFunction* wfun, int narg, char** args);

#endif
