/****
 * @date Created on 2021-03-29 at 16:21:05 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the eigenstate type containing the parameters to compute and store an eigenstate of the wave function
 * in the random Lorentz gas. By eigenstate, we also include resonant state. All the eigenstates of the multiple scattering system
 * are defined by the eigenvalue equation M(k).a = 0, where "M(k)" is the multiple scattering matrix, i.e., inverse transition operator,
 * at wave number "k", and "a" is the eigenvector associated with a zero eigenvalue, mu=0. To determine a nontrivial solution for "a" to the
 * eigenvalue equation, we need to find "k" such that det(M(k))=0.
 ***/
#ifndef _EIGENSTATE_TYPE_H
#define _EIGENSTATE_TYPE_H
#include "dcomplex_type.h"

/**
 * Data contained in an eigenstate in the complex "k" plane.
 */
typedef struct Eigenstate_s {
	
	dcomplex kroot; //Wave number at the position of the eigenstate in the complex "k" plane. Also complex root of det(M(k))=0 for "k".
	dcomplex mu;    //Eigenvalue of the matrix M(k) corresponding to the eigenstate. In principle, this eigenvalue should be as close to zero as possible to ensure the eigenvector relevance.
	dcomplex* eigv; //Complex amplitudes of the eigenstate on each atom (size "med->n"). The atoms are supposed to be stored in the same order as in the medium.
	
	//Fiting parameters:
	int isfitted;   //Boolean flag which 1 if the eigenstate fitting is successful, 0 otherwise.
	double* center; //Center of the eigenstate in the medium (size "med->d"), obtained from non-linear fitting by an exponential packet (if fitted=1), or by weighted average (if fitted=0).
	double locfac;  //Localization factor assuming the shape |a(x)| = exp(-locfac*|x-c| + ampcst), where "c" denotes the eigenstate center.
	double ampcst;  //Amplitude logarithmic constant term of the exponential packet. If the state is normalized to 1, then "ampcst" should be about 0 in case of localization.
	
} Eigenstate;

#endif
