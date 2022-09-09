/****
 * @date Created on 2021-04-08 at 19:35:25 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the eigenstate utilities.
 ***/
#ifndef _EIGENSTATE_UTIL_H
#define _EIGENSTATE_UTIL_H
#include "eigenstate_type.h"
#include "medium_type.h"

void init_eigenstate(Eigenstate* eigst, Medium* med);
void del_eigenstate(Eigenstate* eigst);
void fit_eigenstate_exp(Eigenstate* eigst, Medium* med, int maxit, int verb);
void plot_state_radial(Eigenstate* eigst, Medium* med, int i, FILE* tikzfp);
void plot_state_proj(Eigenstate* eigst, Medium* med, int i, FILE* tikzfp);
void save_eigenstate(Eigenstate* eigst, Medium* med, int i, FILE* fp);
void parse_eigenstate(Eigenstate* eigst, Medium* med, int i, int narg, char** args);

#endif
