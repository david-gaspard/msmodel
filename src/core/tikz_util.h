/****
 * @date Created on 2021-02-18 at 15:46:12 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing utilities to generate TikZ/PGFPlots figures from data.
 ***/
#ifndef _TIKZ_UTIL_H
#define _TIKZ_UTIL_H
#include "curve_type.h"

void init_curve(Curve* crv, int nx, const char* opt);
void del_curve(Curve* crv);
void tikz_plot(FILE* fp, Curve* crv, int nc, char* axisopt, char* precmd);

#endif
