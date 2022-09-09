/****
 * @date Created on 2021-12-21 at 17:12:39 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities to create and manipulate incident wave function structures.
 ***/
#include <stdlib.h>               /* Standard Library for Memory Allocation */
#include <stdio.h>                /* Standard Library for Input and Output */
#include <string.h>               /* Standard Library for String Manipulation */
#include "common_util.h"          /* Import the General Purpose Functions */
#include "real_vector_util.h"     /* Import the Complex Data Parsing functions */
#include "green_bessel.h"         /* Import the Green Functions */
#include "medium_util.h"          /* Import the Medium Utilities */
#include "incident_wave_type.h"   /* Import the Incident Wave Type */

/**
 * Returns the type of incident wave as a string.
 */
char* get_type_iwave(IncidentWave* iwave) {//Private
	switch (iwave->source) {
		case plane:
			return "plane";
		case spherical:
			return "spherical";
		default:
			return "none";
	}
}

/**
 * Sets the type of the given incident wave "iwave" according to the string "str".
 * Returns 1 if the string is not recognized.
 */
int string_to_source(const char* str, IncidentWaveSource* source) {//Private
	if (strcmp(str, "plane") == 0) {
		*source = plane;
	}
	else if (strcmp(str, "spherical") == 0) {
		*source = spherical;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Initialize the incident wave using the given arguments. Returns no failure code 0.
 */
int init_iwave(IncidentWave* iwave, IncidentWaveSource source, dcomplex k) {
	iwave->k = k;
	iwave->source = source;
	sprintf(iwave->desc_latex, "%s $k=(%g%+g\\,{\\rm i})\\,\\varsigma^{-1}$", get_type_iwave(iwave), creal(k), cimag(k));
	sprintf(iwave->desc_ascii, "%s k=%g%+gi", get_type_iwave(iwave), creal(k), cimag(k));
	return 0;
}

/**
 * Copy the content of "iwave_src" to "iwave_dst".
 */
void copy_iwave(IncidentWave* iwave_src, IncidentWave* iwave_dst) {
	init_iwave(iwave_dst, iwave_src->source, iwave_src->k);
}

/**
 * Computes the value of the incident wave at the given position.
 */
dcomplex eval_incident_wave(IncidentWave* iwave, int d, double* pos) {
	switch (iwave->source) {
		case plane:
			return cexp(I*iwave->k*pos[0]); //Incident plane wave in the direction 1x.
		case spherical:
			return free_green(d, iwave->k, norm(d, pos)); //Spherical wave centered at the origin.
		default:
			return 0.; //Default value.
	}
}

/**
 * Computes the value of the total wave function (incident+scattered) at the given position "pos".
 */
dcomplex eval_total_wave(IncidentWave* iwave, Medium* med, dcomplex* a, double* pos) {
	dcomplex sca = 0.;
	double r;
	int i, d = med->d;
	for (i = 0; i < med->n; i++) {//Loop on the scatterers.
		r = distance(d, med->pos+i*d, pos);  //Distance between current point and scatterer position.
		sca += a[i]*free_green(d, iwave->k, r);
	}
	return sca + eval_incident_wave(iwave, d, pos);
}

/**
 * Saves the given incident wave "iwave" to the given file.
 */
void save_iwave(IncidentWave* iwave, FILE* fp) {
	fprintf(fp, "k=%.16g%+.16gi\n", creal(iwave->k), cimag(iwave->k));
	fprintf(fp, "source=%s\n", get_type_iwave(iwave));
}

/**
 * Parse the given incident wave using the arguments stored in "args".
 * The "iwave" pointer is supposed to be already allocated and can be dereferenced.
 */
int parse_iwave(IncidentWave* iwave, int narg, char** args) {//Public
	int fail = 0; //The parsing does not fail by default..
	int nargmin = 2; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments, found %d but expected at least %d...\n", narg, nargmin);
		fail = 1;
	}
	char* value;
	double rek, imk;
	IncidentWaveSource source;
	value = get_value(narg, args, "k");
	if (sscanf(value, "%lg%lgi", &rek, &imk) != 2) {
		printf("[ERROR] Invalid complex wave number '%s'...\n", value);
		fail = 1;
	}
	dcomplex k = rek + I*imk;
	value = get_value(narg, args, "source");
	if (string_to_source(value, &source)) {
		printf("[ERROR] Invalid source '%s'...\n", value);
		fail = 1;
	}
	fail |= init_iwave(iwave, source, k);
	return fail;
}
