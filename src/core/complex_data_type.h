/****
 * @date Created on 2021-02-17 at 11:56:31 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the general complex data structure which is inherited by several other structures (complex_map/complex_cut/fast_cut).
 * The idea is to exploit emulated polymorphism to perform the same computation with all these structures (see ms_core.c).
 ***/
#ifndef _COMPLEX_DATA_TYPE_H
#define _COMPLEX_DATA_TYPE_H
#include "dcomplex_type.h"

typedef struct ComplexData_s {
	int ntot;          //Total number of samples contained in the data, whatever the associated structure (rectangle/triplet/line).
	int nseed;         //Number of seeds combined/averaged together in the data.
	dcomplex* data;    //Array containing "ntot" complex-valued samples in all.
	char title[80];    //Short title of the data set, such as 'Trace-log function' or 'Resonance density'.
	double realtime;   //Total computation time taken to process the data (in seconds).
	dcomplex (*get_arg)(void* self, int i);  //Function returning the complex position in wave number "k" corresponding to the given index "i" in the data (provided at initialization).
} ComplexData; //Suffix "_t" is deprecated, since conflicting with POSIX conventions. Here using capitals instead.

#endif
