/****
 * @date Created on 2021-03-21 at 12:13:05 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities of the color data which gives a pixel color for any given real value between 0 and 1.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "pixel_type.h"
#include "color_data_type.h"
#include "colormaps.h"

/**
 * Initializes the given color data with the name of the color map.
 */
void init_colordata(ColorData* coldata, const char* name) {
	if (strcmp(name, "magma") == 0) {
		coldata->colors = magma_colors;
		coldata->len = magma_colors_len;
	}
	else if (strcmp(name, "hesperia") == 0) {
		coldata->colors = hesperia_colors;
		coldata->len = hesperia_colors_len;
	}
	else if (strcmp(name, "viridis") == 0) {
		coldata->colors = viridis_colors;
		coldata->len = viridis_colors_len;
	}
	else if (strcmp(name, "parula") == 0) {
		coldata->colors = parula_colors;
		coldata->len = parula_colors_len;
	}
	else if (strcmp(name, "temperature") == 0) {
		coldata->colors = temperature_colors;
		coldata->len = temperature_colors_len;
	}
	else {
		fprintf(stderr, "[ERROR] Unknown color scheme '%s', aborting...\n", name);
		exit(EXIT_FAILURE);
	}
	strncpy(coldata->name, name, 49);
}

/**
 * Clamps the given double value to the unit interval [0, 1].
 */
static double clamp(double x) {
	if (x < 0. || isnan(x))
		return 0.;
	else if (x > 1.)
		return 1.;
	return x;
}

/**
 * Assigns the RGB values of the given pixel "pix" according to the real values "r", "g", "b" contained in the interval [0, 1].
 */
void set_pixel(double r, double g, double b, Pixel* pix) {
	pix->red   = (uint8_t)round(255*clamp(r));
	pix->green = (uint8_t)round(255*clamp(g));
	pix->blue  = (uint8_t)round(255*clamp(b));
}

/**
 * Sets the pixel value "pix" (red, green, blue) to the color linearly interpolated from the data
 * in the color data "coldata", and according to the color value "val" in the interval [0, 1].
 */
void pick_color_from_data(ColorData* coldata, double val, Pixel* pix) {
	double r, g, b;  //Real color values in the interval [0, 1].
	double pos = clamp(val)*(coldata->len-1); //Real-valued index position of the color.
	int i = (int)floor(pos); //Find the index of the smallest color of the interpolation interval (index+1 being the largest color).
	if (i == coldata->len-1) {//Check for special case value = 1 for which the smallest color would be index=len-1, i.e., the last element.
		i--;
	}
	pos = pos - i; //The position is now in the interval [0, 1] between two colors.
	r = (1-pos)*coldata->colors[i][0] + pos*coldata->colors[i+1][0];
	g = (1-pos)*coldata->colors[i][1] + pos*coldata->colors[i+1][1];
	b = (1-pos)*coldata->colors[i][2] + pos*coldata->colors[i+1][2];
	set_pixel(r, g, b, pix);
}
