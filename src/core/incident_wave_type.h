/****
 * @date Created on 2021-12-21 at 16:37:30 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the type of incident wave function.
 ***/
#ifndef _INCIDENT_WAVE_TYPE_H
#define _INCIDENT_WAVE_TYPE_H
#include "dcomplex_type.h"
/**
 * Type of incident wave function.
 */
typedef enum IncidentWaveSource_e {plane, spherical} IncidentWaveSource;

/**
 * Structure defining the incident wave.
 */
typedef struct IncidentWave_s {
	
	IncidentWaveSource source;  //Type of incident wave function (typically "plane" or "spherical").
	dcomplex k;                 //Wavenumber used for the incident wave. May be complex in general.
	
	char desc_ascii[80];        //Short description of the incident wave in ASCII format for printing.
	char desc_latex[80];        //Short description of the incident wave in LaTeX format for printing.
	
} IncidentWave;

#endif
