/****
 * @date Created on 2021-02-23 at 18:44:30 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the logarithmic complex cut utilities.
 ***/
#ifndef _LOG_CUT_UTIL_H
#define _LOG_CUT_UTIL_H
#include "log_cut_type.h" /* Import the Complex Cut Type */
#include "medium_type.h"  /* Import the Medium Type */

void init_logcut(LogCut* lcut, int ns, double rek, double imkmin, double imkmax, int nseed, const char* title);
void del_logcut(LogCut* lcut);
void print_param_logcut(LogCut* lcut);
void save_logcut(LogCut* lcut, const char* fname);
void export_tikz_lcut(LogCut* lcut, Medium* med, const char* fname);
void parse_logcut(LogCut* lcut, int narg, char** args);

#endif
