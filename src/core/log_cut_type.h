/****
 * @date Created on 2021-02-23 at 16:57:30 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex cut type with a logarithmic mesh, i.e., exponentially spaced grid points.
 * The interest of this kind of mesh is the ability of exploring several orders of magnitude in the complex k-plane with a few number of grid points.
 * This is especially useful when plotting in log-log scale to identify power-law behaviors on the resonance density.
 ***/
#ifndef _LOG_CUT_TYPE_H
#define _LOG_CUT_TYPE_H
#include "complex_data_type.h"

typedef struct LogCut_s {
	ComplexData;    //Main data of the fast cut implemented by inheritance (must be first in the struct, cannot be a pointer, FMS extension needed).
	double rek;     //Constant real part of the log cut, typically positive and larger than 5.
	double imkmin;  //Start point of the log cut at k[0] = rek + imkmin*I. The point k=rek is the asymptotic point of the logarithmic mesh.
	double imkmax;  //End point of the log cut at k[ntot-1] = rek + imkmax*I. The sign of "imkmax" must be the same of "imkmin".
	dcomplex* k;    //List of the positions where the function should be evaluated. The ordinates are stored in "data", see ComplexData.
} LogCut;

#endif
