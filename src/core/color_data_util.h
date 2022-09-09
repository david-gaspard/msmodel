/****
 * @date Created on 2021-03-21 at 13:29:14 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the color data utilities. This is essentially the function to colorize a pixel from agiven value in the interval [0, 1].
 ***/
#ifndef _COLOR_DATA_UTIL_H
#define _COLOR_DATA_UTIL_H
#include "color_data_type.h"
#include "pixel_type.h"

void init_colordata(ColorData* coldata, const char* name);
void pick_color_from_data(ColorData* coldata, double val, Pixel* pix);

#endif
