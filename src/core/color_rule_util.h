/****
 * @date Created on 2021-03-22 at 18:53:13 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the color rule utilities.
 ***/
#ifndef _COLOR_RULE_UTIL_H
#define _COLOR_RULE_UTIL_H
#include "color_rule_type.h"
#include "pixel_type.h"

void set_color(ColorRule* colrule, double h, Pixel* pix);
char* tikz_colorbar_v1(ColorRule* colrule);
char* tikz_colorbar_v2(ColorRule* colrule);
void save_colorrule(ColorRule* colrule, char* name, FILE* fp);
void parse_colorrule(ColorRule* colrule, char* arg);

#endif
