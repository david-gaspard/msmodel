/***
 * @date Created on 2020-08-08 at 11:34 CEST
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Create and manage the medium consisting of N point-like atoms.
 **/
#include <stdlib.h>                 /* Standard Library for Memory Allocation */
#include <stdio.h>                  /* Standard Library for Input and Output */
#include <string.h>                 /* Standard Library for String Manipulation */
#include <math.h>                   /* Standard Library for Mathematical Functions */
#include "medium_type.h"            /* Import the Medium Data Structure */
#include "common_util.h"            /* Import the General Purpose Functions */
#include "real_vector_util.h"       /* Import the Real Vector Utilities (distance function) */
#include "nru_random.h"             /* Import the Numerical Recipes Random Generator (thread-safe version) */
#include "green_bessel.h"           /* Import the Complex Library and the Green and Bessel functions */
#include "scattering_model_util.h"  /* Import the Scattering Model Utilities */
/**
 * Define some constant macros:
 */
#define PI          3.141592653589793238
#define TWOPI       6.2831853071795864769  //Value of 2*pi.
#define SQRTTWOPI   2.506628274631000502   //Square root of 2*pi.
#define NUMEPS      1.e-16                 //Machine epsilon in double precision.
#define LNMAX       709.7827128933839731   //Maximum value of the argument of a natural exponential function in double precision.

/**
 * Initialization function of an empty medium with the given parameters.
 * This function does not fill the medium with the atoms.
 * The content of the medium must be deleted with the corresponding del_*() function.
 */
void init_medium(Medium* med, int d, int n, ScatteringModel* scmodel, Shape shape, double ratio) {
	if (n < 2) {
		printf("[ERROR] The number of atoms should be strictly greater than 1, aborting...\n");
		exit(EXIT_FAILURE);
	}
	if (n > 5e4) {
		printf("[ERROR] The number of atoms n=%d is too large, aborting...\n", n);
		exit(EXIT_FAILURE);
	}
	if (d < 1) {
		printf("[ERROR] The dimension should be strictly positive, aborting...\n");
		exit(EXIT_FAILURE);
	}
	if (d > 20) {
		printf("[ERROR] The dimension d=%d is too large, aborting...\n", d);
		exit(EXIT_FAILURE);
	}
	if (ratio < 1.) {
		printf("[ERROR] The aspect ratio should be larger or equal to 1, aborting...\n");
		exit(EXIT_FAILURE);
	}
	if (d == 1 && shape == ball) {
		shape = cube;
		printf("[WARN] Changed medium shape 'ball' to 'cube' since 1D.\n");
	}
	med->d = d;
	med->n = n;
	med->ntot = n*d; //Total number of coordinates. Note that "seed" is initialized later when filling the medium.
	med->shape = shape;
	med->ratio = ratio;
	med->pos = (double*)calloc(med->ntot, sizeof(double)); //Use calloc to set the values to zero, just in case all the data are not set when printing a file.
	med->rnd = (Urandom*)calloc(1, sizeof(Urandom));   //Allocate space for a local random generator, for thread safety.
	med->scmodel = (ScatteringModel*)calloc(1, sizeof(ScatteringModel));
	init_scmodel(med->scmodel, scmodel->type, scmodel->np, scmodel->p); //Initializes a deep copy of the given scattering model.
}

/**
 * Initialize a deep copy of the given medium "srcmed" into "destmed".
 */
void init_copy_medium(Medium* srcmed, Medium* dstmed) {
	init_medium(dstmed, srcmed->d, srcmed->n, srcmed->scmodel, srcmed->shape, srcmed->ratio);
}

/**
 * Deletes the content (allocated pointers) within the given medium, not the pointer itself. Do not forget to invoke after use of "init_*" functions.
 */
void del_medium(Medium* med) {
	free(med->scmodel);
	free(med->pos);
	free(med->rnd);
}

/***************************
 * MEDIUM FILLING FUNCTIONS
 ***************************/
/**
 * Fills the given medium "med" with uniformly randomly placed atoms in a d-cube.
 * Ensures the density of atoms is equal to 1, so that each atom occupies a unit volume on average.
 */
void fill_medium_cube(Medium* med) {
	double s = pow(med->n, 1./med->d); //Length of the cube.
	int l;
	for (l = 0; l < med->ntot; l++) {
		med->pos[l] = s*random_double(med->rnd);
	}
}

/**
 * Returns the radius of the medium ensuring unit atomic density as if the medium was a ball. This radius is of course exact for a ball-shaped medium.
 * Since the density is maintained equal to one, this can be used for instance in spherical approximation of the actual medium shape.
 */
double radius_of_ball(Medium* med) {
	return pow(med->n/uball_volume(med->d), 1./med->d);
}

/**
 * Fills the given medium "med" with randomly placed atoms in a d-ball.
 * Ensures the density of atoms is equal to 1, so that each atom occupies a unit volume on average.
 */
void fill_medium_ball(Medium* med) {
	double radius = radius_of_ball(med);  //Radius of the d-ball to ensure unit density.
	double norm2, r2 = radius*radius;
	int i, c, l; //Atom index, coordinate index, total index.
	for (i = 0; i < med->n; i++) {//Loop on the atoms.
		l = i*med->d;
		do {//Rejection method.
			norm2 = 0.0;
			for (c = 0; c < med->d; c++) {
				med->pos[l+c] = radius*(2*random_double(med->rnd)-1);
				norm2 += med->pos[l+c]*med->pos[l+c];
			}
		} while (norm2 > r2);
	}
}

/**
 * Fills the given medium "med" with a regular cubic lattice in arbitrary dimension (integer lattice Z^d).
 * Use a change of base method to fill a cube as equal as possible.
 */
void fill_medium_lattice(Medium* med) {
	int s = (int)ceil(pow(med->n, 1./med->d)); //Integer length of the cube.
	int i, c, val;
	for (i = 0; i < med->n; i++) {
		val = i;
		for (c = 0; c < med->d; c++) {//Change of base method ("val" in base "s").
			med->pos[i*med->d+c] = val%s;
			val /= s;
		}
	}
}

/**
 * Fills the given medium "med" with a gaussian density in arbitrary dimension.
 * The standard deviation of the gaussian density, i.e., the "radius" of the medium,
 * is set such that the density at the center (x=0) is equal to 1 sp^-d by definition of "sp".
 */
void fill_medium_gaussian(Medium* med) {
	double sigma = pow(med->n, 1./med->d)/SQRTTWOPI; //Ensures the density is equal to 1 sp^-d at the center.
	int i;
	for (i = 0; i < med->ntot; i++) {
		med->pos[i] = sigma*random_normal(med->rnd); //Same formula for all coordinates (factorizable distribution).
	}
}

/**
 * Fills the given medium "med" with randomly placed atoms according to the prescriptions given in "medium".
 * Apply a volume-preserving stretching transformation in the first coordinate only if "ratio" is not 1.
 */
void fill_medium(Medium* med, uint64_t seed) {
	init_random(med->rnd, seed);  //Re-initializes the random generator.
	switch (med->shape) {
		case cube:
			fill_medium_cube(med);
			break;
		case ball:
			fill_medium_ball(med);
			break;
		case lattice:
			fill_medium_lattice(med);
			break;
		case gaussian:
			fill_medium_gaussian(med);
			break;
	}
	if (med->ratio != 1. && med->d != 1) {
		double fac2 = pow(med->ratio, -1./med->d); //Scale factor in the (d-1) squeezed directions, f2 = (1/r)^(1/d) < 1.
		double fac1 = med->ratio*fac2; //Scale factor in the stretched direction, f1 = r*(1/r)^(1/d) > 1.
		int l;
		for (l = 0; l < med->ntot; l++) {
			if (l%med->d == 0) {//Stretching transformation in the first coordinate.
				med->pos[l] *= fac1;
			}
			else {
				med->pos[l] *= fac2;
			}
		}
	}
}

/**
 * Converts the given string to the corresponding shape.
 * Returns 1 if the string is not recognized. 
 */
int string_to_shape(const char* str, Shape* shape) {
	if (strcmp(str, "cube") == 0) {
		*shape = cube;
	}
	else if (strcmp(str, "ball") == 0) {
		*shape = ball;
	}
	else if (strcmp(str, "lattice") == 0) {
		*shape = lattice;
	}
	else if (strcmp(str, "gaussian") == 0) {
		*shape = gaussian;
	}
	else {
		return 1;
	}
	return 0;
}

/**
 * Converts the given shape to the corresponding string.
 * Returns the string "none" if the shape is not recognized. 
 */
char* shape_to_string(Shape shape) {
	switch (shape) {
		case cube:
			return "cube";
		case ball:
			return "ball";
		case lattice:
			return "lattice";
		case gaussian:
			return "gaussian";
		default:
			return "none";
	}
}

/**
 * Computes the Euclidean distance between the atoms of index "i" and "j" in the given medium.
 * The indices i, j must be in {0, ..., N-1}. This function does not perform bounds checking.
 */
double atom_distance(Medium* med, int i, int j) {
	return distance(med->d, med->pos + i*med->d, med->pos + j*med->d);
}

/**
 * Returns the maximum distance between two atoms in the medium.
 * Useful to estimate the range of validity of the potential method: kimax=ln(eps)/max_distance(med), with eps=1e-16.
 * This function assumes that the medium is already filled with the correct seed.
 */
double max_distance(Medium* med) {
	double r, rmax = atom_distance(med, 0, 1);
	int i, j;
	for (j = 2; j < med->n; j++) {//Loop over the columns, hence column-major ordering.
		for (i = 0; i < j; i++) {//Loop over the rows in the upper triangle (matrix is symmetric).
			r = atom_distance(med, i, j);
			if (rmax < r)
				rmax = r;
		}
	}
	return rmax;
}

/**
 * Returns th value of "kimax". This is the maximum imaginary part of the wave number "k" allowed in the computations of the MS matrix.
 * It is defined as the value of "k" for which the relative difference between two MS matrix elements exceed 10^16. Note that "kimax" is generally negative.
 * Use the given "seed" to construct the medium and estimate "kimax".
 */
double kimax_value(Medium* med, uint64_t seed) {
	fill_medium(med, seed); //Initialize the medium with the given seed.
	return log(NUMEPS)/max_distance(med);
}

/**
 * Computes the mass center of the given medium "med", as given by the arithmetic average of the atomic positions.
 * The mass center position is written to the array "center" of size "d".
 */
void center_medium(Medium* med, double* center) {
	int i, c;
	for (c = 0; c < med->d; c++) {//Loop on the coordinates of the mass center.
		center[c] = 0.;
		for (i = 0; i < med->n; i++) {//Loop on the atoms.
			center[c] += med->pos[i*med->d + c];
		}
		center[c] /= med->n;
	}
}

/***********************
 * SCATTERING FUNCTIONS
 **********************/
/**
 * Returns the inverse of the single-atom scattering amplitude according to the prescribed model "scmodel".
 * If the scattering model is unknown, then returns a default zero value (unphysical infinite cross section).
 * This function enters the diagonal elements of the multi-scattering matrix.
 */
dcomplex invf(Medium* med, dcomplex k) {
	return invf_scmodel(med->scmodel, med->d, k);
}

/**
 * Returns the single-atom cross section in units of sp^(d-1).
 * Continues analytically wherever k is complex valued.
 */
dcomplex cross_section(Medium* med, dcomplex k) {
	return (conj(1./invf(med, conj(k))) - 1./invf(med, k))/(2*I*k);
}

/**
 * Returns the maximum single-atom cross section allowable at the given wave number.
 * Continues analytically wherever k is complex valued.
 */
dcomplex max_cross_section(Medium* med, dcomplex k) {
	return 8*PI/(k*uball_surface(med->d)*cpow(k/TWOPI, med->d-2));
}

/**
 * Returns the mean free path of the medium at the given real wave number "k", i.e., l=1/(n*sigma).
 * Since the atomic density, "n", is unity (choice of units), this is just the inverse of the cross section.
 */
double mean_free_path(Medium* med, double k) {
	return cabs(1./cross_section(med, k));
}

/**************************
 * INPUT/OUTPUT FUNCTIONS
 *************************/
/**
 * Prints a short summary of the parameters used in medium for the standard output.
 */
void print_param_medium(Medium* med) {
	printf("[INFO] Medium %dD of shape '%s' (ratio=%g) with %d atoms and model %s %s.\n",
		med->d, shape_to_string(med->shape), med->ratio, med->n, med->scmodel->name, med->scmodel->param_ascii);
	//if (verbose) {
	//	int nsugg = (int)round(med->len1)*pow(round(med->len2), med->d-1);
	//	printf("[INFO] Suggested size %d %d, leading to %d atoms with ratio=%g.\n", (int)round(med->len1), (int)round(med->len2), nsugg, round(med->len1)/round(med->len2));
	//}
}

/**
 * Prints the information about the single-atom cross section, especially comparing with the maximum value.
 * Also compares the mean free path and the approximate system radius.
 */
void print_cross_section(Medium* med, dcomplex k) {
	dcomplex cs = cross_section(med, k);
	printf("[INFO] At k = %g%+gi sp^-1, cross section = %g%+gi sp^%d", creal(k), cimag(k), creal(cs), cimag(cs), med->d-1);
	if (cimag(k) == 0.) {//Prints the maximum cross section only for purely real wave number "k".
		dcomplex maxcs = max_cross_section(med, k);
		printf(", maximum = %g%+gi sp^%d (%g%%)", creal(maxcs), cimag(maxcs), med->d-1, 100*cabs(cs/maxcs));
	}
	printf(".\n");
	double radius = radius_of_ball(med);  //Radius of the medium as if it is a d-ball.
	double mfp = 1/cabs(cs);  //Mean free path for a unit density medium.
	printf("[INFO] Mean free path l=%g sp, radius R=%g sp, ratio R/l=%g.\n", mfp, radius, radius/mfp);
}

/**
 * Detects if the given bound min[Im(k)] < 0 on the imaginary part of the wave number is likely to produce roundoff errors while computing on the multi-scattering matrix.
 * Roundoff errors are expected when Im(k) is too negative, or if the maximum atom-atom distance is too large, e.g., in strongly elongated media.
 * Only prints a status message if a possible problem is detected.
 */
void detect_roundoff_error(Medium* med, double imk, int verbose) {
	double kimax = kimax_value(med, 1);
	//double len2 = pow(med->n/med->ratio, 1./med->d); //Approximate length of the medium in the (d-1) small direction(s). L2 = (N/r)^(1/d).
	//double len1 = len2*med->ratio; //Approximate length of the medium in the elongated direction. L1 = r*(N/r)^(1/d)
	//double kimax = -LNMAX/sqrt(len1*len1 + (med->d - 1)*len2*len2);  //Critical minimum value of Im(k) allowed.
	if (imk < kimax) {
		printf("[WARN] Imaginary part %g is too negative. Risk of roundoff errors in the MS matrix. Make it larger than %g.\n", imk, kimax);
	}
	else if (verbose) {
		printf("[INFO] No risk of roundoff errors with Im(k) = %g. The critical region is Im(k) < %g.\n", imk, kimax);
	}
}

/**
 * Exports the medium "med" to a TikZ picture file. In dimensions greater than 3, projects in the first three dimensions, omitting the last coordinates.
 * When the given "seed" is non-zero, then shows an outlook of the medium for this particular seed.
 */
void export_tikz_medium(Medium* med, const char* fname, uint64_t seed) {
	char tikzfname[strlen(fname)+6];
	sprintf(tikzfname, "%s.tikz", fname);
	FILE* tikzfp = fopen(tikzfname, "a");
	fprintf(tikzfp, "\\begin{center}%%\nMedium %dD of shape ``%s'' (ratio $%g$) with $N=%d$ atoms,\\\\ model %s %s\n\\end{center}%%\n",
		med->d, shape_to_string(med->shape), med->ratio, med->n, med->scmodel->name, med->scmodel->param_latex);
	if (seed != 0) {//If nonzero seed, then draws the medium with the given seed.
		fill_medium(med, seed);  //First ensure that the medium is filled with default seed 1.
		fprintf(tikzfp, "\\par\\begin{tikzpicture}%%\n\\begin{axis}[%%\n");
		fprintf(tikzfp, "\ttitle={Seed $%lu$},\n", seed);
		fprintf(tikzfp, "\twidth=0.5\\textwidth,\n\tenlargelimits=true,\n"); //\taxis equal image,\n
		int i; //Atom index.
		if (med->d == 1) {//Dimension 1.
			fprintf(tikzfp, "\txlabel={$x~(\\varsigma)$},\n");
			fprintf(tikzfp, "]\n\\addplot[black, mark=*, only marks, mark options={black, scale=0.5}, ycomb, black!40] coordinates {%%\n");
			for (i = 0; i < med->n; i++) {//Loop on atoms.
				fprintf(tikzfp, "(%g,1) ", med->pos[i]);
			}
		}
		else if (med->d == 2) {//Dimension 2.
			fprintf(tikzfp, "\txlabel={$x~(\\varsigma)$},\n\tylabel={$y~(\\varsigma)$},\n");
			fprintf(tikzfp, "]\n\\addplot[black, mark=*, only marks, mark options={black, scale=0.5}] coordinates {%%\n");
			for (i = 0; i < med->n; i++) {//Loop on atoms.
				fprintf(tikzfp, "(%g,%g) ", med->pos[2*i], med->pos[2*i+1]);
			}
		}
		else {//Dimension 3 and more.
			fprintf(tikzfp, "\tview={15}{30},\n");
			fprintf(tikzfp, "\txlabel={$x~(\\varsigma)$},\n\tylabel={$y~(\\varsigma)$},\n\tzlabel={$z~(\\varsigma)$},\n");
			fprintf(tikzfp, "]\n\\addplot3[mark=*, only marks, mark options={black, scale=0.5}, ycomb, black!40] coordinates {%%\n");
			for (i = 0; i < med->n; i++) {//Loop on atoms. Ignore coordinates of index larger than 3.
				fprintf(tikzfp, "(%g,%g,%g) ", med->pos[med->d*i], med->pos[med->d*i+1], med->pos[med->d*i+2]);
			}
		}
		fprintf(tikzfp, "\n};\n\\end{axis}\n\\end{tikzpicture}\n");
	}
	fclose(tikzfp);
}

/**
 * Saves the medium data, excluding the atomic positions, to the file "fname".
 */
void save_medium(Medium* med, const char* fname) {
	FILE* fp = fopen(fname, "a");
	fprintf(fp, "[medium]\n");
	fprintf(fp, "dimension=%d\n", med->d);
	fprintf(fp, "natom=%d\n", med->n);
	save_scmodel(med->scmodel, fp);
	fprintf(fp, "shape=%s\n", shape_to_string(med->shape));
	fprintf(fp, "ratio=%g\n\n", med->ratio);
	fclose(fp);
}

/**
 * Saves the positions of the atoms stored in "med" in CSV format (comma-separated value) to a file of automatic name.
 * Overwrite any existing file of the same name. To be used for testing purpose only.
 */
void save_medium_points(Medium* med) {
	char fname[50];
	sprintf(fname, "medium_%s_seed%lu.tab", shape_to_string(med->shape), med->rnd->seed);
	FILE* fp = fopen(fname, "w"); //Overwrite mode.
	int i, c;
	for (i = 0; i < med->n; i++) {//Loop on the atoms.
		for (c = 0; c < med->d; c++) {//Loop on the coordinates of each atom.
			fprintf(fp, "%.16g", med->pos[i*med->d + c]);
			if (c < med->d-1) {
				fprintf(fp, ","); //No space.
			}
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	printf("[INFO] Medium points saved to file '%s'...\n", fname);
}

/**
 * Parses the medium item "med" with the arguments of "args".
 * The "med" pointer is supposed to be already allocated and can be dereferenced.
 */
void parse_medium(Medium* med, int narg, char** args) {
	int nargmin = 5; //Minimum number of argument.
	if (narg < nargmin) {
		printf("[ERROR] Too few arguments in medium, found %d but expected at least %d, aborting...\n", narg, nargmin);
		exit(EXIT_FAILURE);
	}
	int d, n;      //Mandatory fields.
	double ratio;  //Scattering length.
	Shape shape;   //Medium shape.
	char* value = get_value(narg, args, "dimension");
	if (sscanf(value, "%d", &d) != 1) {
		printf("[ERROR] Invalid dimension '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "natom");
	if (sscanf(value, "%d", &n) != 1) {
		printf("[ERROR] Invalid number of atoms '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "shape");
	if (string_to_shape(value, &shape)) {
		printf("[ERROR] Invalid shape '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "ratio");
	if (sscanf(value, "%lg", &ratio) != 1) {
		printf("[ERROR] Invalid aspect ratio '%s', aborting...\n", value);
		exit(EXIT_FAILURE);
	}
	value = get_value(narg, args, "model");
	if (!value[0]) {//Model not found
		printf("[ERROR] Missing scattering model, aborting...\n");
		exit(EXIT_FAILURE);
	}
	ScatteringModel scmodel;
	parse_scmodel(&scmodel, value);
	init_medium(med, d, n, &scmodel, shape, ratio);  //Call a constructor-like function init_medium(Medium* med, int d, int n, ...).
}
