/****
 * @date Created on 2021-09-16 at 10:12:49 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the scattering model utilities.
 ***/
#ifndef _SCATTERING_MODEL_UTIL_H
#define _SCATTERING_MODEL_UTIL_H
#include "dcomplex_type.h"
#include "scattering_model_type.h"  /* Import the Scattering Model Structure */

int init_scmodel(ScatteringModel* scmodel, ScatteringType type, int np, double* p);
void del_scmodel(ScatteringModel* scmodel);
dcomplex invf_scmodel(ScatteringModel* scmodel, int d, dcomplex k);
void save_scmodel(ScatteringModel* scmodel, FILE* fp);
void parse_scmodel(ScatteringModel* scmodel, char* arg);

#endif
