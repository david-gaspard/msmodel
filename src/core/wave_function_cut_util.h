/****
 * @date Created on 2021-12-21 at 20:33:34 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the wave function cut utilities.
 ***/
#ifndef _WAVE_FUNCTION_CUT_UTIL_H
#define _WAVE_FUNCTION_CUT_UTIL_H
#include "wave_function_cut_type.h"   /* Import the Wave Function Cut Type */
#include "medium_type.h"              /* Import the Medium Type */

void del_wfcut(WaveFunctionCut* wfcut);
void print_param_wfcut(WaveFunctionCut* wfcut);
void get_pos_wfcut(WaveFunctionCut* wfcut, int l, double* pos);
void store_data_wfcut(WaveFunctionCut* wfcut, double* repsi_samples, double* impsi_samples, double* density_samples);
void export_tikz_wfcut(WaveFunctionCut* wfcut, Medium* med, const char* fname);
void save_wfcut(WaveFunctionCut* wfcut, const char* fname);
void parse_wfcut(WaveFunctionCut* wfcut, int narg, char** args);

#endif
