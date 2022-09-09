/****
 * @date Created on 2021-09-14 at 15:31:07 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex histogram type describing the complex eigenvalues of the multiple scattering matrix.
 ***/
#ifndef _COMPLEX_HISTOGRAM_TYPE_H
#define _COMPLEX_HISTOGRAM_TYPE_H
#include "dcomplex_type.h"
#include "domain_type.h"
#include "color_rule_type.h"
/**
 * Type of the multiple scattering matrix for which the eigenvalue histogram is being computed.
 * "msmatrix": Usual MS matrix defined as M(k) = F(k)^-1 - G(k).
 * "normalized": Dimensionless MS matrix defined as N(k) = i - G(k)/I(k,0), corresponding to phase shift kept equal to delta=90°.
 */
typedef enum MatrixType_e {msmatrix, normalized} MatrixType;
/**
 * Complex histogram type for the matrix eigenvalues
 */
typedef struct ComplexHistogram_s {
	
	MatrixType mtype;    //Type of the matrix being computed (either "msmatrix" or "normalized").
	int xbin;            //Horizontal number of bins.
	int ybin;            //Vertical number of bins.
	int ntot;            //Total number of bins (nx*ny).
	Domain* dom;         //Rectangular domain boundaries. Center coordinates of the extreme pixels.
	int dom_ok;          //Flag indicating if the domain ranges xmin:xmax and ymin:ymax are set up (0=not set up, 1=already set up).
	
	dcomplex k;          //Complex wave number at which the eigenvalues histogram of the MS matrix are computed.
	int nseed;           //Number of different atomic geometries, or seeds, that should be averaged/combined in the mu-plane histogram.
	dcomplex avgz;       //Average of the complex numbers entered in the histogram.
	dcomplex sigmaz;     //Standard deviation in real and imaginary parts around the average.
	double* data;        //Histogram of the eigenvalue distribution of the MS matrix at given wave number "k". Number of eigenvalues in each bin. Natural row-major ordering.
	
	ColorRule* colrule;  //Color scheme used for image output.
	double realtime;     //Total computation time (in seconds).
	char title[80];      //Title of the complex histogram (optional).
	
	//Marginal distribution of the imaginary part of the eigenvalues:
	int mbin;      //Number of bins used for the marginal distribution of the imaginary part.
	double mmin;   //Minimum abscissa of the histogram. Lower bound of the smallest bin. Computed internally.
	double mmax;   //Maximum abscissa of the histogram. Upper bound of the largest bin. Computed internally.
	double* mdata; //Histogram data. Number of imaginary part of eigenvalues in each bin (size "nbin").
	int mrange_ok; //Flag indicating if the range of the marginal distribution is set up.
	
} Chistogram;

#endif
