/****
 * @date Created on 2021-02-17 at 12:26:59 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing a few utilities for the general complex data structure.
 ***/
#ifndef _COMPLEX_DATA_UTIL_H
#define _COMPLEX_DATA_UTIL_H
#include "complex_data_type.h"

void save_cdata(ComplexData* cdata, FILE* fp);
void parse_cdata(ComplexData* cdata, char* arg);

#endif
