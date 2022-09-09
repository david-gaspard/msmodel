/****
 * @date Created on 2021-05-17 at 18:26:00 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to plot several complex functions in the complex "k" plane.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "pixel_type.h"
#include "common_util.h"
#include "complex_vector_util.h"
#include "medium_util.h"
#include "complex_map_util.h"
#include "ms_matrix_util.h"

#define PI         3.1415926535897932384626  //Defines the Pi constant.
#define TWOPI      6.2831853071795864769253  //Defines the 2*Pi constant.
#define PId3       1.0471975511965977461542  //Defines the Pi/3 constant.

#define RAYNUM     32
#define MINVAL     0.5

/**
 * Converts a given HSV color ([-pi,pi], [0,1], [0,1]) to a pixel defined as a standard RGB color ([0,255], [0,255], [0,255]).
 */
void hsv_to_rgb(double hsv[3], Pixel* pix) {
	double hue = hsv[0]/PId3;  //Gets the color hue and maps it to interval [-3,3].
	double sat = hsv[1];       //Gets the color saturation.
	double val = hsv[2]*255;   //Gets the color value and multiplies by 255 so that it runs over [0,255].
	double r, g, b;
	switch ((int)floor(hue)) {  //Convert color hue to RGB:
	case -3:  //hue in [-3, -2]
		r = 0.0;
		g = -hue - 2.0;
		b = 1.0;
		break;
	case -2:  //hue in [-2, -1]
		r = hue + 2.0;
		g = 0.0;
		b = 1.0;
		break;
	case -1:  //hue in [-1, 0]
		r = 1.0;
		g = 0.0;
		b = -hue;
		break;
	case 0:  //hue in [0, 1]
		r = 1.0;
		g = hue;
		b = 0.0;
		break;
	case 1:  //hue in [1, 2]
		r = -hue + 2.0;
		g = 1.0;
		b = 0.0;
		break;
	default:  //hue in [2, 3]
		r = 0.0;
		g = 1.0;
		b = hue - 2.0;
		break;
	}
	//Accounts for the saturation and value:
	pix->red   = (uint8_t)(val*(1.0-sat + sat*r));
	pix->green = (uint8_t)(val*(1.0-sat + sat*g));
	pix->blue  = (uint8_t)(val*(1.0-sat + sat*b));
}

/**
 * Colorize the given pixel "pix" from the value of the given complex number "z".
 * Use a version of the domain color code developped by E. Wegert.
 */
void colorize(dcomplex z, Pixel* pix) {
	double hsv[3];         //Double HSV color defined in (hue=[-pi,pi], sat=[0,1], val=[0,1]);
	double m = cabs(z);    //Compute the absolute value of z (modulus).
	if (m < INFINITY && m == m) {   //Ensure the modulus to be neither Infinity nor NaN.
		hsv[0] = carg(z);                //Compute the color hue in the interval [-pi,pi] from the argument of z.
		//double amod = RAYNUM*(hsv[0] + PI)/TWOPI;
		//amod = amod - floor(amod);
		double amod = 1.;
		
		hsv[1] = 1.0;       //Always maximum saturation in this color map.
		double lmod = RAYNUM*log(m)/TWOPI;
		lmod = lmod - floor(lmod);                     //Compute the fractional part.
		hsv[2] = MINVAL + (1.-MINVAL)*(amod + lmod)/2; //Prevent the saturation value from reaching a totally white color (sat=0).
	}
	else {//Values for Infinity and NaN.
		hsv[0] = 0.;
		hsv[1] = 0.;
		hsv[2] = 1.;
	}
	hsv_to_rgb(hsv, pix);
}

/**
 * Export the given complex map in PPM format, assuming a complex domain coloring.
 */
void export_ppm_cmap(ComplexMap* cmap, FILE* fp) {
	fprintf(fp, "P6\n%d %d\n255\n", cmap->nx, cmap->ny);  //Header with 'P6' for binary PPM file format.
	Pixel pix;
	int l;
	for (l = 0; l < cmap->ntot; l++) {//Loop on pixels to write.
		colorize(cmap->data[l], &pix);
		fwrite(&pix, 1, sizeof(Pixel), fp);
	}
}

/**
 * Fill the complex map with the values computed from the function (serial computation).
 */
void compute_cmap(ComplexData* cdata, dcomplex (*f)(dcomplex)) {
	double begin = omp_get_wtime();  //Captures the beginning time (omp-shared constant).
	char msg[] = "Computing function";
	dcomplex z;
	int l, dl = cdata->ntot/2000;
	if (dl == 0) dl = 1;
	for (l = 0; l < cdata->ntot; l++) {//Loop on the points to compute.
		z = cdata->get_arg(cdata, l);
		cdata->data[l] = f(z);
		if (l%dl == 1)
			print_progress(l, cdata->ntot, omp_get_wtime()-begin, msg);
	}
	cdata->realtime = omp_get_wtime()-begin;  //Saves the total computation time.
	print_end(l, cdata->realtime, msg);
}

/**
 * Main function.
 */
int main(int argc, char** argv) {
	
	//Initialize a medium:
	char* med_args[] = {
		"dimension=3",
		"natom=40",
		"model=hardsphere",
		"alpha=0.1",
		"shape=ball",
		"ratio=1.0"
	};
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	parse_medium(med, sizeof(med_args)/sizeof(med_args[0]), med_args);
	uint64_t seed = 1; //Seed of the medium.
	fill_medium(med, seed);
	print_param_medium(med);
	
	//Initialize a complex map:
	char* cmap_args[] = {
		"title=Test Function",
		"nseed=1", //Not used here...
		"xrange=5:8",
		"yrange=-2:0",
		//"xrange=-3:3",
		//"yrange=-3:3",
		"xsample=600",
		"color=temperature div 5"  //Color will be changed anyway...
	};
	ComplexMap* cmap = (ComplexMap*)calloc(1, sizeof(ComplexMap));
	parse_complexmap(cmap, sizeof(cmap_args)/sizeof(cmap_args[0]), cmap_args);
	
	dcomplex* matrix = (dcomplex*)calloc(med->n*med->n, sizeof(dcomplex));
	dcomplex* vec = (dcomplex*)calloc(med->n, sizeof(dcomplex));
	constant_unit_cvector(med->n, vec); //Reset the vector to normalized vector.
	
	dcomplex mineigval_1(dcomplex k) {//Minimum eigenvalue function.
		int ivit_maxit = 50;  //Maximum number of iterations of the inverse power iteration.
		double ivit_toler = 1e-10;  //Tolerance of the relative distance between two successive vectors obtained by inverse iterations.
		int ivit_verb = 0;  //Verbosity level of the inverse iteration function (0=quiet, 1=verbose).
		dcomplex mu;
		build_ms_matrix(med, k, matrix); //Builds the MS matrix at "kroot" in prevision of the inverse iteration.
		inverse_iteration(med->n, matrix, vec, &mu, ivit_maxit, ivit_toler, ivit_verb); //Execute a few inverse iterations to compute the eigenvector/eigenvalue.
		return mu;
	};
	
	dcomplex mineigval_2(dcomplex k) {//Minimum eigenvalue function.
		build_ms_matrix(med, k, matrix); //Builds the MS matrix at "kroot" in prevision of the inverse iteration.
		eigvals(med->n, matrix, vec);
		dcomplex mu = vec[0];
		int i;
		for (i = 1; i < med->n; i++) {//Search for the minimum eigenvalue.
			if (cabs(vec[i]) < cabs(mu)) {
				mu = vec[i];
			}
		}
		return mu;
	};
	
	//Defines the function that will be plotted:
	dcomplex func(dcomplex k) {
		//return k; //Identity to show the bare color code.
		//return dk_log_det_ms(med, k); //Logarithmic derivative of det(M(k)).
		build_ms_matrix(med, k, matrix);
		//return cexp(med->n*tracelog_lu(med->n, matrix)); //Determinant det(M(k)) itself, for small "n" only.
		//return mineigval_1(k);
		//return mineigval_2(k);
		return 1./traceinv(med->n, matrix);
	};
	
	//Fill the complex map with the values computed from the function:
	compute_cmap((ComplexData*)cmap, func);
	
	//Draw the complex map to some PPM file.
	FILE* fp = fopen("untitled.ppm", "w");
	export_ppm_cmap(cmap, fp);
	fclose(fp);
	
	//Frees memory:
	del_medium(med);
	del_complexmap(cmap);
	free(med);
	free(cmap);
	free(vec);
	free(matrix);
	return 0;
}
