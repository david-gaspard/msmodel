/****
 * @date Created on 2021-07-05 at 18:08:22 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the random function type containing the samples of a random variable which depends on a certain abscissa.
 * This abscissa is typically the wave number "k", and the random variable is typically the total cross section for a given incident direction, or the density of states.
 ***/
#ifndef _RANDOM_FUNCTION_TYPE_H
#define _RANDOM_FUNCTION_TYPE_H
/**
 * Type of random function to be computed (crsec=cross_section, dos=density_of_states).
 */
typedef enum RandomFunctionType_e {crsec, dos} RandomFunctionType;
/**
 * Structure containing the output data.
 */
typedef struct RandomFunction_s {
	RandomFunctionType type;  //Type of random function to be computed (either cross section or density of states).
	int nbin;                 //Number of separated bins from "xmin" to "xmax" (included).
	int nseed;                //Total number of atomic geometries (seeds) used to make the statistics.
	int ntot;                 //Total size of the data array (nbin*ns).
	double xmin;              //Minimum abscissa, corresponding to the bin of index 0.
	double xmax;              //Maximum abscissa, corresponding to the bin of index nbin-1.
	double* data;             //Array containing all the samples of the function (size ntot = nbin*nseed).
	                          //The first level is for the bins, and the second for the samples within the bin.
	char title[80];           //Title of the random function (optional).
	double realtime;          //Real computation time (in seconds).
} RandomFunction;

#endif
