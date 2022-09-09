/****
 * @date Created on 2020-11-30 at 14:42:22 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex cut type. This represents a cartesian linear mesh of (n+2) samples long and 3 samples wide.
 ***/
#ifndef _COMPLEX_CUT_TYPE_H
#define _COMPLEX_CUT_TYPE_H
#include "complex_data_type.h"

typedef struct ComplexCut_s {
	ComplexData;          //Main data of the fast cut implemented by inheritance (must be first in the struct, cannot be a pointer, FMS extension needed).
	int ns;               //Number of effective samples along the cut. All the effective samples are surrounded by eight neighbors.
	//int ntot;             //Total number of samples, ntot = 3*(ns+2).
	dcomplex zmin;        //Start point of the mesh. 
	dcomplex zmax;        //End point of the mesh.
	dcomplex dz;          //Complex step between consecutive samples, dz = (zmax-zmin)/(ns-1).
	double h;             //Length of the step between consecutive samples, h = |dz| = |zmax-zmin|/(ns-1).
	//int nseed;            //Number of seeds combined/averaged together in the complex cut.
	//char title[80];       //Short title of the complex cut, such as 'Trace-log function' or 'Resonance density'.
	//double realtime;      //Total computation time taken to fill the complex cut (in seconds).
	//dcomplex* data;       //List of the samples in triplet-major ordering (running over the individual triplets from left to right in reading order, then from zmin to zmax).
} ComplexCut;

#endif
