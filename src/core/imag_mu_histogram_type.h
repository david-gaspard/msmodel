/****
 * @date Created on 2021-07-31 at 13:42:47 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the histogram of imaginary parts of the complex eigenvalues "mu" of the multiple scattering matrix.
 * This histogram uses a logarithmic scale.
 ***/
#ifndef _IMAG_MU_HISTOGRAM_TYPE_H
#define _IMAG_MU_HISTOGRAM_TYPE_H

typedef struct ImagMuHistogram_s {
	
	int nseed;         //Number of different atomic geometries.
	double k;          //Wave number at which the histogram is computed.
	int nbin;          //Number of bins of the histogram data.
	double xmin;       //Minimum abscissa of the histogram. Lower bound of the smallest bin. Computed internally.
	double xmax;       //Maximum abscissa of the histogram. Upper bound of the largest bin. Computed internally.
	double* data;      //Histogram data. Number of imaginary part of eigenvalues in each bin (size "nbin").
	double realtime;   //Total computation time in seconds.
	char title[80];    //Short title of the job.
	
	double avgx;       //True average of the imaginary parts of "mu". Computed internally.
	double varx_true;  //True variance of the imaginary part of "mu", computed directly from the samples.
	double varx_expc;  //Expected variance of the imaginary parts of "mu", computed from Tr[G⁺G]/(4N).
	
} ImagMuHistogram;

#endif
