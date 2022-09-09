/****
 * @date Created on 2021-09-15 at 19:43:35 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the utilities of the scattering model, especially the functions to compute the inverse scattering amplitude.
 ***/
#include <stdlib.h>                 /* Standard Library for Memory Allocation */
#include <stdio.h>                  /* Standard Library for Input and Output */
#include <string.h>                 /* Standard Library for String Manipulation */
#include <math.h>                   /* Standard Library for Mathematical Functions */
//#include "dcomplex_type.h"          /* Import the Complex Type */
#include "common_util.h"            /* Import the General Purpose Functions */
#include "green_bessel.h"           /* Import the Complex Library and the Green and Bessel functions */
#include "real_vector_util.h"       /* Import the Real Vector Utilities */
#include "scattering_model_type.h"  /* Import the Scattering Model Structure */
/**
 * Define some constant macros:
 */
#define PI          3.141592653589793238   //Famous pi constant.
#define TWOPI       6.2831853071795864769  //Value of 2*pi.

/**
 * Returns the expected number of parameters for the given type of scattering model.
 */
int expected_nparam(ScatteringType type) {
	switch (type) {
		case hardsphere:
			return 1;
		case softsphere:
			return 1;
		case resonant:
			return 2;
		case maximum:
			return 0;
		default:
			return 0;
	}
}

/**
 * Initialization function of a scattering model with the given parameters.
 * The content of the medium must be deleted with the corresponding del_*() function.
 */
void init_scmodel(ScatteringModel* scmodel, ScatteringType type, int np, double* p) {
	int np_expc = expected_nparam(type);
	if (np != np_expc) {
		printf("[ERROR] Unexpected number of scattering parameters. Given %d but expected %d, aborting...\n", np, np_expc);
		exit(EXIT_FAILURE);
	}
	switch (type) {
		case hardsphere:
			sprintf(scmodel->dir, "a%g", p[0]);  //Short name of the model used for the directory. Example: "a0.1"
			sprintf(scmodel->name, "hardsphere"); //Name of the model for ASCII/LaTeX without parameter. Example: "hardsphere"
			sprintf(scmodel->param_ascii, "alpha=%gsp", p[0]); //Parameters of the model for ASCII only. Example: "alpha=0.1sp"
			sprintf(scmodel->param_latex, "$\\alpha=%g\\,\\varsigma$", p[0]); //Parameters of the model for LaTeX only. Example: "$\\alpha=0.1\\,\\varsigma$"
			break;
		case softsphere:
			sprintf(scmodel->dir, "b%g", p[0]);
			sprintf(scmodel->name, "softsphere");
			sprintf(scmodel->param_ascii, "alpha=%gsp", p[0]);
			sprintf(scmodel->param_latex, "$\\alpha=%g\\,\\varsigma$", p[0]);
			break;
		case resonant:
			sprintf(scmodel->dir, "r%g", p[0]);
			sprintf(scmodel->name, "resonant");
			sprintf(scmodel->param_ascii, "p=%g%+gi sp^-1", p[0], p[1]);
			sprintf(scmodel->param_latex, "$p=(%g%+g\\,{\\rm i})\\varsigma^{-1}$", p[0], p[1]);
			break;
		case maximum:
			sprintf(scmodel->dir, "max");
			sprintf(scmodel->name, "maximum");
			scmodel->param_ascii[0] = '\0';
			scmodel->param_latex[0] = '\0';
			break;
	}
	scmodel->type = type;
	scmodel->np = np;
	scmodel->p = (double*)calloc(np, sizeof(double));
	int i;
	for (i = 0; i < np; i++) {//Copy the parameters.
		scmodel->p[i] = p[i];
	}
}

/**
 * Deletes the content (allocated pointers) within the given scattering model, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_scmodel(ScatteringModel* scmodel) {
	free(scmodel->p);
}

/**
 * Returns the inverse of the single-atom scattering amplitude in our custom soft-sphere model.
 * NB: This model is superseded by the hard-sphere model, which has the same energy behavior at low energy, but is more physical.
 */
dcomplex invf_ss(ScatteringModel* scmodel, int d, dcomplex k) {//Private function
	double alpha = scmodel->p[0];
	if (d == 2) {
		//return I/4 - clog(alpha*k)/TWOPI; //Old convention, but then k=i/alpha is not the single-atom pole of F(k).
		return clog(I/(alpha*k))/TWOPI;
	}
	int drec = d;
	double prefac, frec = TWOPI*alpha*alpha;
	if (d%2 == 0) {
		prefac = 1./4;
	}
	else {
		prefac = alpha/2;
	}
	while (drec > 2) {
		drec -= 2;
		prefac /= frec*drec;
	}
	//printf("[TEST] Prefac(d=%d, alpha=%f) = %f.\n", d, alpha, prefac);
	//Note that k=i/alpha must be the single-atom pole of F(k). TODO: Find convenient generalization for d=1 and d>3.
	return prefac*(1. + I*cpow(alpha*k, d-2));
}

/**
 * Returns the inverse of the single-atom scattering amplitude corresponding to hard spheres of radius "alpha". This model has no bound states.
 * The conservation of probability at real "k" is guaranteed, even for large values of "alpha", because the scatterers are treated as points.
 */
dcomplex invf_hs(ScatteringModel* scmodel, int d, dcomplex k) {//Private function
	double alpha = scmodel->p[0];
	dcomplex i0 = uball_surface(d)*cpow(k/TWOPI, d-2)/(8*PI);
	return -i0*free_green(d, k, alpha)/free_green_imag(d, k, alpha);
}

/**
 * Returns the inverse of the single-atom scattering amplitude corresponding to the resonant sphere model (pure Breit-Wigner profile).
 * The resonance pole has the exact position p = pre + I*pim, where "pim" should be negative.
 */
dcomplex invf_res(ScatteringModel* scmodel, int d, dcomplex k) {//Private function
	double pre = scmodel->p[0]; //Position of the single-atom resonance in the "k" plane (with negative imaginary part), at which F(k) is infinite.
	double pim = scmodel->p[1];
	dcomplex i0 = uball_surface(d)*cpow(k/TWOPI, d-2)/(8*PI);
	return i0*(I - (k - pre)/pim);
}

/**
 * Returns the inverse of the single-atom scattering amplitude corresponding to the maximum model.
 */
dcomplex invf_max(int d, dcomplex k) {//Private function
	dcomplex i0 = uball_surface(d)*cpow(k/TWOPI, d-2)/(8*PI);
	return I*i0;
}

/**
 * Returns the inverse of the single-atom scattering amplitude for the prescribed model.
 */
dcomplex invf_scmodel(ScatteringModel* scmodel, int d, dcomplex k) {//Public function
	switch (scmodel->type) {
	case hardsphere:
		return invf_hs(scmodel, d, k);
	case softsphere:
		return invf_ss(scmodel, d, k);
	case resonant:
		return invf_res(scmodel, d, k);
	case maximum:
		return invf_max(d, k);
	default:
		return 0.;
	}
}

/**
 * Save the given scattering model to the given file.
 */
void save_scmodel(ScatteringModel* scmodel, FILE* fp) {
	fprintf(fp, "model=%s", scmodel->name);
	int i;
	for (i = 0; i < scmodel->np; i++) {
		fprintf(fp, " %.16g", scmodel->p[i]);
	}
	fprintf(fp, "\n");
}

/**
 * Converts the given string "str" to the corresponding scattering model type.
 * Returns 1 if the string is not recognized. 
 */
int string_to_sctype(const char* str, ScatteringType* type) {
	if (strcmp(str, "hardsphere") == 0) {
		*type = hardsphere;
	}
	else if (strcmp(str, "softsphere") == 0) {
		*type = softsphere;
	}
	else if (strcmp(str, "resonant") == 0) {
		*type = resonant;
	}
	else if (strcmp(str, "maximum") == 0) {
		*type = maximum;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Parses the given scattering model with the argument "arg".
 * The "scmodel" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_scmodel(ScatteringModel* scmodel, char* arg) {
	ScatteringType type;
	char str[150];
	int ic;  //Current number of parsed characters. Used as an index.
	sscanf(arg, "%148s %n", str, &ic);
	if (string_to_sctype(str, &type)) {
		printf("[ERROR] Invalid scattering model '%s', aborting...\n", str);
		exit(EXIT_FAILURE);
	}
	int np = count_real_data(arg+ic);
	double* p = (double*)calloc(np, sizeof(double));
	parse_real_data(np, p, arg+ic);
	init_scmodel(scmodel, type, np, p);
	free(p);
}
