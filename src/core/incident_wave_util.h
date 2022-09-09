/****
 * @date Created on 2021-12-21 at 17:30:26 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the incident wave utilities.
 ***/
#ifndef _INCIDENT_WAVE_UTIL_H
#define _INCIDENT_WAVE_UTIL_H
#include "incident_wave_type.h"   /* Import the Incident Wave Type */
#include "medium_type.h"          /* Import the Medium Type */

void copy_iwave(IncidentWave* iwave_src, IncidentWave* iwave_dst);
dcomplex eval_incident_wave(IncidentWave* iwave, int d, double* pos);
dcomplex eval_total_wave(IncidentWave* iwave, Medium* med, dcomplex* a, double* pos);
void save_iwave(IncidentWave* iwave, FILE* fp);
int parse_iwave(IncidentWave* iwave, int narg, char** args);

#endif
