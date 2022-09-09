/****
 * @date Created on 2021-07-10 at 12:27:42 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the type containning the differential cross section data.
 ***/
#ifndef _DIFF_CROSS_SECTION_TYPE_H
#define _DIFF_CROSS_SECTION_TYPE_H

typedef struct DiffCrossSection_s {
	int nseed;          //Number of independent seeds used to generate the random medium. Each atomic geometric leads to a whole curve of differential cross section.
	uint64_t iseed;     //Index of the first seed to be used (included). Useful to change the atomic geometry.
	double k;           //Wave number "k" at which the differential cross section is computed (units of sp^-1).
	double thmin;       //Minimum angle in degree (typically 0), corresponding to the bin of index 0.
	double thmax;       //Maximum angle in degree (typically 180), corresponding to the bin of index nbin-1.
	int nbin;           //Number of separated bins in the angular interval thmin:thmax (including bounds).
	int nq;             //Number of quantiles to represent the statistical distribution of the differential cross section in each bin (typically 5: 0/4, 1/4, 2/4, 3/4, 4/4).
	int ntot;           //Total size of the data array (nbin*nq). Used only in internals.
	double* data;       //Array containing the quantiles of the function (size ntot = nbin*nq).
	                    //The first level is for the bins, and the second for the quantiles within the bin.
	double* mean;       //Array containing the mean curve as a single number per bin (size nbin).
	int nc;             //Number of desired sample curves for the first seeds.
	double* curve;       //Array containing the curves obtained with the "nc" first seeds (size nc*nbin). The first level is for the bin, and the second for the curve.
	char title[80];     //Title of the random function (optional).
	double realtime;    //Real computation time (in seconds).
} DiffCrossSection;

#endif
