/****
 * @date Created on 2020-08-12 at 11:14:43 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex map type.
 ***/
#ifndef _COMPLEX_MAP_TYPE_H
#define _COMPLEX_MAP_TYPE_H
#include "complex_data_type.h"
#include "color_rule_type.h"

typedef struct ComplexMap_s {
	ComplexData;         //Main data of the complex map implemented by inheritance (must be first in the struct, cannot be a pointer, FMS extension needed).
	int nx;              //Number of samples in the horizontal direction.
	int ny;              //Number of samples in the vertical direction.
	double xmin;         //Minimum x value (included), left border.
	double xmax;         //Maximum x value (included), right border.
	double ymin;         //Minimum y value (included), lower border.
	double ymax;         //Maximum y value (included), upper border.
	double dx;           //Horizontal step between consecutive samples, dx = (xmax-xmin)/(nx-1).
	double dy;           //Vertical step between consecutive samples, dy = (ymax-ymin)/(ny-1). Generally equal to dx.
	ColorRule* colrule;  //Color scheme used for image output.
} ComplexMap;

#endif
