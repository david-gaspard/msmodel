/****
 * @date Created on 2021-04-05 at 18:06:10 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the eigenstate list type which defines a simulation to search for a many eigenstates.
 ***/
#ifndef _EIGENSTATE_LIST_TYPE_H
#define _EIGENSTATE_LIST_TYPE_H
#include "domain_type.h"
#include "medium_type.h"
#include "eigenstate_type.h"
/**
 * Enumerates the methods that can be used to find the eigenstates.
 */
typedef enum EigenstateMethod_e {determinant, mineigval, invtraceinv} EigenstateMethod;

/**
 * Structure containing "nstate" eigenstates along with the root-finding settings.
 */
typedef struct EigenstateList_s {
	
	Medium* med;   //Reference to the medium which contains the dimension, number of atoms, and atomic positions.
	uint64_t seed; //Seed of the random medium for which the eigenstates will be computed.
	
	//Root-finding parameters:
	EigenstateMethod method; //Numerical method used to find the eigenstates in the complex "k" plane.
	int ntarget;   //User-defined number of targets exploited in the root-finding process to search for the eigenstates in the domain "dom".
	Domain* dom;   //Complex domain of research for the eigenstates in the "k" plane (xmin, xmax, ymin, ymax).
	int maxit;     //Maximum number of secant root-finding steps in the complex "k" plane.
	double toler;  //Tolerance on the relative error between consecutive steps in the complex "k" plane.
	int verb;      //Verbosity level of the root finder (0=quiet, 1=verbose, 2=debug).
	
	//Resulting data:
	int nstate;           //Final number of states found by the root-finding process. This is used to allocate space for "eigtab". We have always: nstate <= ntarget.
	int nfitted;          //Total number of successfully fitted eigenstates.
	Eigenstate** eigtab;  //Data associated to the different eigenvectors (size "nstate"). More flexible to have "Eigenstate** eigtab" in order to swap eigenstates by reference more easily (needed for "uniq" and "sort").
	double realtime;      //Total computation time taken to find all the roots and the eigenstates (in seconds).
	
} EigenstateList;

#endif
