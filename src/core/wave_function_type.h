/****
 * @date Created on 2021-11-27 at 18:19:14 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the wave function structure.
 ***/
#ifndef _WAVE_FUNCTION_TYPE_H
#define _WAVE_FUNCTION_TYPE_H
#include <stdint.h>  /* Standard Library for Fixed-Size Integers */
#include "domain_type.h"
#include "incident_wave_type.h"
#include "color_rule_type.h"

/**
 * Structure of the wave function object.
 */
typedef struct WaveFunction_s {
	
	int nseed;                  //Number of seeds used when averaging over the random configurations of the scatterer positions.
	uint64_t iseed;             //Initial seed used (included).
	IncidentWave* iwave;        //Incident wave function, containing the type (typically "plane" or "spherical") and the wave number.
	
	int nx;                     //Number of samples in the horizontal direction.
	int ny;                     //Number of samples in the vertical direction.
	Domain* dom;                //Domain where the wave function is defined.
	double dx;                  //Horizontal step between consecutive samples, dx = (xmax-xmin)/(nx-1).
	double dy;                  //Vertical step between consecutive samples, dy = (ymax-ymin)/(ny-1). Generally equal to dx.
	
	double* function;           //2D array with the real part of the wave function, averaged over the random configurations (size nx*ny, row-major format).
	double* density;            //2D array with the square modulus of the wave function, averaged over the random configurations (size nx*ny, row-major format).
	
	ColorRule* colrule1;        //Color scheme used for image output.
	ColorRule* colrule2;        //Color scheme used for image output.
	double realtime;            //Total computation time.

} WaveFunction;

#endif
