/****
 * @date Created on 2021-02-20 at 11:51:53 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the curve type to be used for PGFPlots/TikZ plotting purposes.
 ***/
#ifndef _CURVE_TYPE_H
#define _CURVE_TYPE_H

typedef struct Curve_s {
	int nx;        //Number of points.
	double* x;     //Abscissas of points.
	double* y;     //Ordinates of points.
	char opt[80];  //TikZ plotting options.
} Curve;

#endif
