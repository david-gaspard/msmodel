/****
 * @date Created on 2021-04-01 at 16:10:33 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the utilities to manipulate eigenstate structures.
 ***/
#ifndef _EIGENSTATE_LIST_UTIL_H
#define _EIGENSTATE_LIST_UTIL_H
#include "eigenstate_list_type.h"
#include "medium_type.h"

void alloc_eigenstate_list(EigenstateList* eigls, int nstate);
void del_eigenstate_list(EigenstateList* eigls);
char* eigmethod_to_string(EigenstateMethod method);
void check_eigenstate_list(EigenstateList* eigls);
void export_tikz_eigenstate_list(EigenstateList* eigls, const char* fname);
void save_eigenstate_list(EigenstateList* eigls, const char* fname);
void parse_eigenstate_list(EigenstateList* eigls, Medium* med, int narg, char** args);

#endif
