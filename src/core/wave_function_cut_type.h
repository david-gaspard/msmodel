/****
 * @date Created on 2021-12-21 at 16:33:51 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the wave function cut type.
 ***/
#ifndef _WAVE_FUNCTION_CUT_TYPE_H
#define _WAVE_FUNCTION_CUT_TYPE_H
#include <stdint.h>  /* Standard Library for Fixed-Size Integers */
#include "incident_wave_type.h"

/**
 * Structure of the wave function object.
 */
typedef struct WaveFunctionCut_s {
	
	int nseed;             //Number of seeds used when averaging over the random configurations of the scatterer positions.
	uint64_t iseed;        //Initial seed used (included).
	IncidentWave* iwave;   //Incident wave function, containing the type (typically "plane" or "spherical") and the wave number.
	
	int nbin;              //Number of samples in the horizontal direction.
	double xmin;           //Lower bound of the interval of coordinates.
	double xmax;           //Lower bound of the interval of coordinates.
	double dx;             //Horizontal step between consecutive samples, dx = (xmax-xmin)/(nx-1).
	
	int ntot;              //Total number of double values in the following array. Typically, ntot = 3*nbin.
	double* repsi;         //1D array with the real part of the wave function, format [mean, quantile_1/4, quantile_3/4] (size 3*nbin).
	double* impsi;         //1D array with the imaginary part of the wave function, format [mean, quantile_1/4, quantile_3/4] (size 3*nbin).
	double* density;       //1D array with the square modulus of the wave function, format [mean, quantile_1/4, quantile_3/4] (size 3*nbin).
	
	double realtime;       //Total computation time.

} WaveFunctionCut;

#endif



