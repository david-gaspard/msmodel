/****
 * @date Created on 2021-09-15 at 19:33:00 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the type of the scattering model for the single-atom collision.
 ***/
#ifndef _SCATTERING_MODEL_H
#define _SCATTERING_MODEL_H
/**
 * Define the scattering model of the point atoms.
 */
typedef enum ScatteringType_e {hardsphere, softsphere, resonant, maximum} ScatteringType;

/**
 * Scattering model structure describing the single-atom collision.
 */
typedef struct ScatteringModel_s {
	
	ScatteringType type;  //Type of scattering model (see enum above).
	
	int np;     //Number of parameters. of the scattering model (only depends on the type).
	double* p;  //List of real parameters of the scattering model.
	
	char dir[50];  //Short name of the model used for the directory. Example: "a0.1"
	char name[150];         //Name of the model for ASCII/LaTeX without parameter. Used for I/O. Example: "hardsphere"
	char param_ascii[150];  //Parameters of the model for ASCII only. Example: "alpha=0.1sp"
	char param_latex[150];  //Parameters of the model for LaTeX only. Example: "$\\alpha=0.1\\,\\varsigma$"
	
} ScatteringModel;

#endif
