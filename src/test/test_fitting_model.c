/****
 * @date Created on 2021-02-21 at 16:44:31 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the test functions of the fitting model.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include "fitting_model_util.h"  /* Imports the Fitting Model Algorithms */

/**
 * Fills the given array "x" with "nx" regularly spaced values from "xmin" to "xmax".
 * The array "x" is assumed to have enough space.
 */
void fill_range(int nx, double xmin, double xmax, double* x) {
	double h = (xmax - xmin)/(nx - 1);
	int i;
	for (i = 0; i < nx; i++) {
		x[i] = xmin + h*i;
	}
}

/**
 * Test the Gauss-Newton algorithm with a straight line.
 */
void test_linear() {
	printf("====== Testing Damped Gauss-Newton method with straight line ======\n");
	FittingModel* fm = (FittingModel*)calloc(1, sizeof(FittingModel));
	int id = 1, np = 2, nx = 10;  //Polynomial least square with two parameters of 10 grid points.
	double xmin = 0, xmax = 1, x[nx];
	fill_range(nx, xmin, xmax, x);
	double y[] = {0.190434, 0.766852, 2.17372, 2.81268, 3.84482, 4.73944, 6.02535, 6.83869, 7.9376, 9.12028};
	
	init_model(fm, id, np, nx, x, y);
	fit_model(fm);
	//Expected -0.03923364826204505 + 8.968441571765503*x.
	
	del_model(fm);
	free(fm);
}


/**
 * Test the Gauss-Newton algorithm with a Padé approximant.
 */
void test_pade() {
	printf("====== Testing Damped Gauss-Newton method with Padé approximant ======\n");
	FittingModel* fm = (FittingModel*)calloc(1, sizeof(FittingModel));
	int id = 2, np = 5, nx = 31;  //Padé least square of 31 grid points.
	double xmin = 0, xmax = 10, x[nx];
	fill_range(nx, xmin, xmax, x);
	double y[] = {3.00952, 1.0815, 0.768452, 0.664585, 0.625173, 0.606052, 0.600253, 0.594224, 0.593367, 0.594596, 0.591818, 0.59552, 0.594312, 0.594107, 0.599347, 0.599345, 0.600839, 0.59963, 0.598757, 0.599049, 0.602053, 0.600063, 0.601826, 0.604932, 0.603617, 0.606808, 0.60438, 0.60654, 0.606599, 0.604962, 0.607267};
	
	init_model(fm, id, np, nx, x, y);
	
	//Difficult initial parameters for the algorithm:
	fm->param[0] = 2.0;
	fm->param[1] = 1.0;
	fm->param[2] = 5.0;
	fm->param[3] = 10.0;
	fm->param[4] = 2.0;
	
	fit_model(fm);
	//Expected p0=3.00951, p1=2.08429, p2=4.87417, p3=6.1747, p4=7.78678.
	
	del_model(fm);
	free(fm);
}

/**
 * Test the Gauss-Newton algorithm with a sine wave.
 */
void test_sine_wave() {
	printf("====== Testing Damped Gauss-Newton method with sine wave ======\n");
	FittingModel* fm = (FittingModel*)calloc(1, sizeof(FittingModel));
	int id = 7, np = 3, nx = 31;  //Sine wave least square with 31 grid points.
	double xmin = 0, xmax = 10, x[nx];
	fill_range(nx, xmin, xmax, x);
	double y[] = {0.190434, -0.0252365, 0.580452, 0.400467, 0.587961, 0.605469, 0.976404, 0.833215, 0.932126, 1.07134, 0.693121, 0.892339, 0.541496, 0.255234, 0.494215, 0.195098, 0.0472534, -0.359903, -0.712063, -0.918314, -0.817741, -1.17376, -1.11067, -0.867315, -1.01671, -0.674045, -0.847934, -0.531829, -0.395472, -0.406504, -0.0164609};
	
	init_model(fm, id, np, nx, x, y);
	
	//Difficult initial parameters for the algorithm:
	fm->param[0] = 1.0;
	fm->param[1] = 1.0;
	fm->param[2] = 2.0;
	
	fit_model(fm);
	//Expected A=-0.981664, k=0.630068, phi=3.0836.
	
	del_model(fm);
	free(fm);
}


/**
 * Test the Gauss-Newton algorithm with the piecewise function used to fit
 * the characteristic function, Tr[ln(M)]/N, along a vertical complex cut.
 */
void test_piecewise() {
	printf("====== Testing Damped Gauss-Newton method with piecewise function ======\n");
	FittingModel* fm = (FittingModel*)calloc(1, sizeof(FittingModel));
	int id = 3, np = 4, nx = 31;  //Piecewise least square of 31 grid points.
	double xmin = 0, xmax = 10, x[nx];
	fill_range(nx, xmin, xmax, x);
	double y[] = {0.0132006, -0.000918788, 0.0126433, 0.000608882, 0.00168004, -0.00183241, 0.00769775, 0.0139919, 0.0504242, 0.105397, 0.15693, 0.239023, 0.312163, 0.395279, 0.502432, 0.597065, 0.701196, 0.800293, 0.904098, 1.015, 1.13782, 1.24667, 1.37037, 1.50069, 1.61837, 1.7528, 1.87036, 2.00464, 2.13345, 2.25803, 2.39687};
	
	init_model(fm, id, np, nx, x, y);
	fit_model(fm);
	//Expected x0=1.96555, a0=-0.310637, a1=0.4967, b1=0.98229.
	
	del_model(fm);
	free(fm);
}

/**
 * Main function for testing purpose.
 */
int main(int argc, char** argv) {
	test_linear();
	test_pade();
	test_sine_wave();
	test_piecewise();
	return 0;
}
