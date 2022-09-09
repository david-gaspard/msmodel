/****
 * @date Created on 2021-03-21 at 11:59:18 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the color rule type used to colorize the complex maps.
 ***/
#ifndef _COLOR_RULE_TYPE_H
#define _COLOR_RULE_TYPE_H
#include "color_data_type.h"

/**
 * Defines the type of color rule (sequential=from dark to light, diverging=white at zero height).
 */
typedef enum ColorType_e {sequential, diverging} ColorType;

/**
 * Defines a color rule which contains all the parameters to colorize any height value "h" between "hmin" and "hmax".
 */
typedef struct ColorRule_s {
	ColorData* data;  //Gives a pixel color for each real between 0 and 1.
	ColorType type;   //Type of color rule (either "diverging" or "sequential").
	double hmin;      //Minimum height value of the plotted function.
	double hmax;      //Maximum height value of the plotted function.
	double ctr;       //Contrast parameter whose interpretation depends on the color type.
	int rev;          //Color reverse flag, 0=no reverse, 1=reverse.
	double (*colorvalue)(struct ColorRule_s* self, double h);  //Private method to compute the color value associated with the given height.
} ColorRule;

#endif
