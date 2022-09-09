/****
 * @date Created on 2020-09-09 at 12:02:21 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the configuration type, i.e., the whole framework of the program.
 ***/
#ifndef _CONFIG_TYPE_H
#define _CONFIG_TYPE_H
#include "medium_type.h"
#include "complex_map_type.h"
#include "complex_cut_type.h"
#include "log_cut_type.h"
#include "chistogram_type.h"
#include "eigenstate_list_type.h"
#include "random_function_type.h"
#include "diff_cross_section_type.h"
#include "imag_mu_histogram_type.h"
#include "wave_function_type.h"
#include "wave_function_cut_type.h"

typedef struct Config_s {
	int nmed;    //Number of mediums, must be 1.
	int ncmap;   //Number of complex maps (k plane).
	int nccut;   //Number of complex cuts (k plane).
	int nlcut;   //Number of logarithmic cuts (k plane).
	int nchist;  //Number of complex histograms (mu eigenvalues).
	int neigls;  //Number of eigenstates list objects.
	int nrfun;   //Number of random functions (real k semi-axis).
	int ndcs;    //Number of differential cross sections (scattering angles).
	int nihist;  //Number of histograms of Im(mu) (mu eigenvalues).
	int nwfun;   //Number of wave functions.
	int nwfcut;  //Number of wave function cuts.
	Medium* med;
	ComplexMap* cmap;
	ComplexCut* ccut;
	LogCut* lcut;
	Chistogram* chist;
	EigenstateList* eigls;
	RandomFunction* rfun;
	DiffCrossSection* dcs;
	ImagMuHistogram* ihist;
	WaveFunction* wfun;
	WaveFunctionCut* wfcut;
} Config;

#endif
