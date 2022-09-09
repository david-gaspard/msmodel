/****
 * @date Created on 2020-08-08 at 19:22 CEST
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Provides the medium structure.
 ***/
#ifndef _MEDIUM_TYPE_H
#define _MEDIUM_TYPE_H
#include "scattering_model_type.h"  //Import the scattering model type.
#include "nru_random_type.h"        //Import the random generator internal state.
/**
 * Define the shape of the medium.
 */
typedef enum Shape_e {cube, ball, lattice, gaussian} Shape;

/**
 * Define the medium structure :
 */
typedef struct Medium_s {//Contains all the geometrical information about the medium.
	int d;              //Dimension of space, i.e., number of coordinates per atom.
	int n;              //Total number of atoms, i.e., point-like scattering centers. Also the volume of the medium since each atom occupies a volume of 1 (sp^d=1).
	int ntot;           //Total number of coordinates, i.e., just the product n*d.
	Shape shape;        //Shape of the medium (see enum here above).
	double ratio;       //Aspect ratio of the medium in one direction (the first coordinate). Typically, ratio > 1 to get an elongated system.
	
	ScatteringModel* scmodel;  //Scattering model used for the point collisions with the atoms. This structure also contains the scattering parameters.
	
	double alpha;       //Scattering length "alpha" of the atoms in units of the mean interatomic spacing (sp=1).
	Urandom* rnd;       //Store the internal state of the random generator for thread safety.
	double* pos;        //List of all the coordinates of the atoms (length = n*d). Assumed form: (atom1_x atom1_y atom1_z ...) (atom2_x atom2_y atom2_z ...) (atom3_x atom3_y atom3_z ...) ...
} Medium;

#endif
