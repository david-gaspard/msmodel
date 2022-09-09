/****
 * @date Created on 2021-02-17 at 12:08:30 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing a few utilities for the general complex data structure.
 ***/
#include <stdlib.h>    /* Standard Library */
#include <stdio.h>     /* Standard Library for Input and Output */
#include <string.h>    /* Standard Library for String Manipulation */
#include "complex_vector_util.h"
#include "complex_data_type.h"

/**
 * Saves the complex data stored in "cdata" to the file "fp".
 */
void save_cdata(ComplexData* cdata, FILE* fp) {
	fprintf(fp, "data=");
	int l;
	for (l = 0; l < cdata->ntot; l++) {
		fprintf(fp, "%.16lg%+.16lgi ", creal(cdata->data[l]), cimag(cdata->data[l]));
	}
	fprintf(fp, "\n");
}

/**
 * Parses the first complex numbers from the given data string "arg" to the array cdata->data.
 * Assumes that cdata->data is already allocated with enough place. This somewhat depends on other initializations.
 * If the argument "arg" is empty, then then function does nothing.
 */
void parse_cdata(ComplexData* cdata, char* arg) {
	if (arg != NULL && arg[0] != '\0') {
		parse_complex_data(cdata->ntot, cdata->data, arg);
	}
}

//?? Write additional functions to simplify the code:
//?? init_cdata(ComplexData* cdata, int ntot, int nseed, const char* title);
//?? del_cdata(ComplexData* cdata);
//?? export_tikz_cdata(ComplexData* cdata, FILE* fp, const char* opt); //Insert \addplot[opt] ...; to the given TikZ file.
