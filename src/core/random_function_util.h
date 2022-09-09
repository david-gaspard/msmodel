/****
 * @date Created on 2021-07-06 at 17:39:53 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 ***/
#ifndef _RANDOM_FUNCTION_UTIL_H
#define _RANDOM_FUNCTION_UTIL_H
#include "random_function_type.h"
#include "medium_type.h"

void del_random_function(RandomFunction* rfun);
char* random_function_type_to_string(RandomFunctionType type);
void print_param_random_function(RandomFunction* rfun);
double get_arg_rfun(RandomFunction* rfun, int i);
void sort_random_function(RandomFunction* rfun);
void export_tikz_random_function(RandomFunction* rfun, Medium* med, const char* fname);
void save_random_function(RandomFunction* rfun, const char* fname);
void parse_random_function(RandomFunction* rfun, int narg, char** args);

#endif
