/****
 * @date Created on 2021-03-21 at 19:08:04 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the color rule utilities.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "color_rule_type.h"
#include "color_data_util.h"

/**
 * Returns the color value in the interval [0, 1] associated with the given height "h" for a "sequential" colormap.
 * As a reminder, a sequential colormap goes from dark colors for color value 0, up to light colors for values near 1.
 * This function assumes that "hmin" and "hmax" are already initialized in "colrule".
 */
double colorvalue_seq(ColorRule* colrule, double h) {
	double val, x = (h - colrule->hmin)/(colrule->hmax - colrule->hmin);  //Normalize the height to the interval [0, 1].
	if (fabs(colrule->ctr) < 1e-16) {//Approximate for very small "ctr".
		val = x + x*(1-x)*colrule->ctr/2;
	}
	else {
		val = log1p(expm1(colrule->ctr)*x)/colrule->ctr;
	}
	return colrule->rev ? 1-val : val;
}

/**
 * Returns the color value in the interval [0, 1] associated with the given height "h" for a "diverging" colormap.
 * As a reminder, a diverging colormap passes to a special color (such as white) when the height is zero (h=0).
 * This function assumes that "hmin" and "hmax" are already initialized in "colrule".
 */
double colorvalue_div(ColorRule* colrule, double h) {
	double val, x = h/fmax(fabs(colrule->hmin), fabs(colrule->hmax));  //Normalize the height to the interval [-1, 1] without shifting in order to preserve zero.
	if (fabs(colrule->ctr) < 1e-8) {//Approximate for small "ctr".
		val = x + x*(1 - x*x)*colrule->ctr*colrule->ctr/6;
	}
	else {
		val = asinh(sinh(colrule->ctr)*x)/colrule->ctr;
	}
	val = (val + 1)/2;  //Convert from [-1, 1] to [0, 1].
	return colrule->rev ? 1-val : val;
}

/**
 * Initializes the color rule. Do not forget to call the corresponding del_*() function to free memory.
 */
void init_colorrule(ColorRule* colrule, const char* name, ColorType type, double ctr, int rev) {
	colrule->data = (ColorData*)calloc(1, sizeof(ColorData));
	init_colordata(colrule->data, name);
	colrule->type = type;
	colrule->ctr = ctr;
	colrule->rev = rev;
	switch (type) {
		case sequential:
			colrule->colorvalue = colorvalue_seq;
			break;
		case diverging:
			colrule->colorvalue = colorvalue_div;
			break;
		default:
			fprintf(stderr, "[ERROR] Unknown color type, aborting...\n");
			exit(EXIT_FAILURE);
	}
}

/**
 * Frees the memory space occupied by the color rule. Call for each init_*().
 * WARN: Strangely, this causes double free error.
 */
//void del_colorrule(ColorRule* colrule) {
//	free(colrule->data);
//}

/**
 * Sets the pixel "pix" to the color associated with the height "h".
 * Accounts for the reverse flag "rev" here.
 */
void set_color(ColorRule* colrule, double h, Pixel* pix) {
	double val = colrule->colorvalue(colrule, h); //Color value.
	pick_color_from_data(colrule->data, val, pix);
}

/**
 * Finds the height value "h" corresponding to the given color value "val" in the interval [0, 1].
 * This function allows the colorbar to be drawn without applying the nonlinear transformation defined by the function colorvalue(height).
 * This function uses the damped secant method to solve for the height value.
 * The function also assumes the bounds "hmin" and "hmax" are correctly initialized in the color rule "colrule".
 * Returns 1 when the root "h" is valid, 0 otherwise.
 */
int height_from_colorvalue(ColorRule* colrule, double val, double* h) {
	const int maxit = 30, maxss = 10;   //Maximum number of iterations, typically 20-30 for Newton-like methods, and maximum sub-steps.
	double toler = 1e-13;  //Tolerance of the relative error in the result.
	double h0 = colrule->hmin, h1 = colrule->hmax, dh, hn, cn;
	double c0 = colrule->colorvalue(colrule, h0) - val;
	double c1 = colrule->colorvalue(colrule, h1) - val;
	int s, ss;
	for (s = 1; s <= maxit; s++) {//Loop on secant iterations.
		dh = c1*(h1 - h0)/(c1 - c0);
		for (ss = 1; ss <= maxss; ss++) {//Perform sub-steps with damped corrections "dh" to ensure convergence.
			hn = h1 - dh;
			cn = colrule->colorvalue(colrule, hn) - val;
			if (fabs(cn) < fabs(c1)) //If the function is smaller, then stop here the sub-steps.
				break;
			dh /= 2;
		}
		h0 = h1; c0 = c1;
		h1 = hn; c1 = cn;
		//printf("#%2d | ss=%2d, h1=%.12g, dh=%g, c1=%g.\n", s, ss, h1, dh, c1);
		if (fabs(dh) < toler*fabs(h1)) {
			//printf("[INFO] Stopped at step %d...\n", s);
			break;
		}
	}
	*h = h1;
	return fabs(c1) < toler; //If the final function value is small enough, then declare the root as valid.
}

/**
 * Returns the TikZ code for the color bar corresponding to the given color rule.
 * The string should be freed externally.
 * @version v1 First version of the color bar with ticks at fixed positions.
 */
char* tikz_colorbar_v1(ColorRule* colrule) {
	double hmin = colrule->hmin, hmax = colrule->hmax;
	if (hmin == hmax) {
		fprintf(stderr, "[ERROR] Invalid minimum and maximum height values (hmin=hmax=%g), aborting...\n", hmin);
		exit(EXIT_FAILURE);
	}
	char* str = (char*)calloc(30*colrule->data->len + 300, sizeof(char));  //Allocate string with enough space.
	int s = 0;  //Current character position.
	s += sprintf(str+s, "\tcolorbar,\n\tcolormap={%s}{", colrule->data->name); //Gives the generic name.
	Pixel pix;
	double h, y, c, cmin = colrule->colorvalue(colrule, hmin), cmax = colrule->colorvalue(colrule, hmax);
	int i, n = colrule->data->len;
	for (i = 0; i < n; i++) {//Loop linearly on several color values to draw the color bar.
		y = ((double)i)/(n - 1);  //Coordinate along the colorbar in [0, 1].
		c = cmin + (cmax - cmin)*y; //Linear color value.
		pick_color_from_data(colrule->data, c, &pix);
		s += sprintf(str+s, "rgb255=(%d,%d,%d) ", pix.red, pix.green, pix.blue);
	}
	s += sprintf(str+s, "},\n\tcolorbar style={%%\n\t\tytick={");
	n = 5;  //Show "n" ticks along the colorbar, for instance [0, 1/4, 1/2, 3/4, 1].
	for (i = 0; i < n; i++) {//Loop on the ticks positions along the colorbar.
		y = ((double)i)/(n - 1);  //Coordinate along the colorbar in [0, 1].
		s += sprintf(str+s, "%g", y);
		if (i != n-1) 
			s += sprintf(str+s, ",");
	}
	if (colrule->type == diverging) {//If diverging colormap, then also shows zero.
		c = colrule->colorvalue(colrule, 0); //Color value associated with zero.
		y = (c - cmin)/(cmax - cmin); //Position of zero on the colorbar.
		s += sprintf(str+s, ",%g", y);
	}
	s += sprintf(str+s, "},\n\t\tyticklabels={");
	for (i = 0; i < n; i++) {//Loop on the tick labels of the colorbar.
		y = ((double)i)/(n - 1);  //Coordinate along the colorbar in [0, 1].
		c = cmin + (cmax - cmin)*y; //Linear color value.
		if (height_from_colorvalue(colrule, c, &h)) {//If the height value corresponding to "c" is found, then print it.
			s += sprintf(str+s, "$%.4g$", h);
		}
		else {//If on the contrary, the height value is not found, then prints a question mark.
			s += sprintf(str+s, "?");
		}
		if (i != n-1) 
			s += sprintf(str+s, ",");
	}
	if (colrule->type == diverging) {//If diverging colormap, then also shows zero.
		s += sprintf(str+s, ",$0$");
	}
	s += sprintf(str+s, "},\n\t},\n");
	return str;
}

/**
 * Returns the TikZ code for the color bar corresponding to the given color rule.
 * The string should be freed externally.
 * @version v2 Second version of the color bar with ticks at powers of ten.
 */
char* tikz_colorbar_v2(ColorRule* colrule) {
	double hmin = colrule->hmin, hmax = colrule->hmax;
	if (hmin == hmax) {
		fprintf(stderr, "[ERROR] Invalid minimum and maximum height values (hmin=hmax=%g), aborting...\n", hmin);
		exit(EXIT_FAILURE);
	}
	char* str = (char*)calloc(30*colrule->data->len + 500, sizeof(char));  //Allocate string with enough space.
	int s = 0;  //Current character position.
	s += sprintf(str+s, "\tcolorbar,\n\tcolormap={%s}{", colrule->data->name); //Gives the generic name.
	Pixel pix;
	double y, c, cmin = colrule->colorvalue(colrule, hmin), cmax = colrule->colorvalue(colrule, hmax);
	int i, j, n = colrule->data->len;
	for (i = 0; i < n; i++) {//Loop linearly on several color values to draw the color bar.
		y = ((double)i)/(n - 1);  //Coordinate along the colorbar in [0, 1].
		c = cmin + (cmax - cmin)*y; //Linear color value.
		pick_color_from_data(colrule->data, c, &pix);
		s += sprintf(str+s, "rgb255=(%d,%d,%d) ", pix.red, pix.green, pix.blue);
	}
	s += sprintf(str+s, "},\n\tcolorbar style={%%\n\t\tytick={");
	double hl[] = {-1e4, -1e3, -100, -10, -1, 0, 1, 10, 100, 1e3, 1e4};
	int nhl = sizeof(hl)/sizeof(hl[0]);
	for (i = 0; i < nhl; i++) {//Loop on the positions of the ticks (powers of ten).
		c = colrule->colorvalue(colrule, hl[i]); //Color value associated with zero.
		y = (c - cmin)/(cmax - cmin); //Position of zero on the colorbar.
		if (0. < y && y < 1.) {//If colobar coordinate lis in the expected interval [0,1].
			s += sprintf(str+s, "%g,", y);
		}
	}
	if (str[s-1] == ',') s--; //Removes the last comma.
	s += sprintf(str+s, "},\n\t\tyticklabels={");
	for (i = 0; i < nhl; i++) {//Loop on the powers of ten.
		c = colrule->colorvalue(colrule, hl[i]); //Color value associated with zero.
		y = (c - cmin)/(cmax - cmin); //Position of zero on the colorbar.
		if (0. < y && y < 1.) {//If colobar coordinate lis in the expected interval [0,1].
			s += sprintf(str+s, "$%g$,", hl[i]);
		}
	}
	if (str[s-1] == ',') s--; //Removes the last comma.
	s += sprintf(str+s, "},\n\t\tminor ytick={");
	for (i = 0; i < nhl; i++) {//Loop on the powers of ten.
		for (j = 2; j < 10; j++) {//Loop on the integer minor ticks.
			if (hl[i] != 0.) {//Zero height is irrelevant because already drawn before.
				c = colrule->colorvalue(colrule, j*hl[i]); //Color value associated with zero.
				y = (c - cmin)/(cmax - cmin); //Position of zero on the colorbar.
				if (0. < y && y < 1.) {//If colobar coordinate lis in the expected interval [0,1].
					s += sprintf(str+s, "%g,", y);
				}
			}
		}
	}
	if (str[s-1] == ',') s--; //Removes the last comma.
	s += sprintf(str+s, "},\n\t},\n");
	return str;
}

/**
 * Converts the given string to the corresponding color type.
 * Returns 1 if the string is not recognized. 
 */
int string_to_colortype(const char* str, ColorType* type) {
	if (strcmp(str, "seq") == 0) {
		*type = sequential;
	}
	else if (strcmp(str, "div") == 0) {
		*type = diverging;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Returns the string corresponding to the given color type
 */
char* colortype_to_string(ColorType type) {
	switch (type) {
		case sequential:
			return "seq";
		case diverging:
			return "div";
		default:
			return "?";
	}
}

/**
 * Saves the given color rule to the file "fp" using the given "name" (typically "color").
 */
void save_colorrule(ColorRule* colrule, char* name, FILE* fp) {
	fprintf(fp, "%s=%s %s %g", name, colrule->data->name, colortype_to_string(colrule->type), colrule->ctr);
	if (colrule->rev) {
		fprintf(fp, " rev");
	}
	fprintf(fp, "\n");
}

/**
 * Parse the color rule "colrule" with the given argument string "arg". Consume the given argument "arg".
 * Expected format: <colorname> <type[seq|div]> <ctr> [rev]
 */
void parse_colorrule(ColorRule* colrule, char* arg) {
	if (arg[0] == '\0') {//If no argument, then use default parameters.
		init_colorrule(colrule, "magma", sequential, 3, 0);
		printf("[INFO] Initializing color scheme with default map '%s'...\n", colrule->data->name);
		return;
	}
	char name[50], typestr[30], revstr[10];
	ColorType type = sequential;  //Default type of colormap.
	double ctr;   //Contrast parameter.
	int rev = 0;  //Default reverse flag value.
	int nc;       //Number of parsed characters.
	sscanf(arg, "%49s %n", name, &nc);
	arg += nc;
	sscanf(arg, "%29s %n", typestr, &nc);
	if (string_to_colortype(typestr, &type)) {
		printf("[ERROR] Invalid type of colormap '%s', aborting...\n", typestr);
		exit(EXIT_FAILURE);
	}
	arg += nc;
	if (sscanf(arg, "%lg %n", &ctr, &nc) != 1) {
		printf("[ERROR] Invalid contrast parameter in '%s', aborting...\n", arg);
		exit(EXIT_FAILURE);
	}
	arg += nc;
	if (sscanf(arg, "%9s", revstr) == 1 && strcmp(revstr, "rev") == 0) {
		rev = 1;
		printf("[INFO] Color reversal enabled.\n");
	}
	init_colorrule(colrule, name, type, ctr, rev);
}
