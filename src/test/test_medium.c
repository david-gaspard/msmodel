/****
 * @date Created on 2020-09-20 at 17:15:51 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing testing functions for the construction of the medium and the generation of atomic positions.
 ***/
#include <stdio.h>
#include <stdlib.h>
#include "medium_util.h"

/**
 * Roughly test the shape of the medium.
 * @deprecated Should be moved to the corresponding test code.
 */
void test_medium_1(Medium* med) {
	double L1_max = med->pos[0];  //Boundaries of the medium.
	double L1_min = L1_max;
	double L2_max = med->pos[1];
	double L2_min = L2_max;
	double C1 = 0.0;  //Center of mass of the medium.
	double C2 = 0.0;
	int l;
	for (l = 0; l < med->ntot; l++) {
		if (l%med->d == 0) {
			if (med->pos[l] > L1_max) {
				L1_max = med->pos[l];
			}	
			if (med->pos[l] < L1_min) {
				L1_min = med->pos[l];
			}
			C1 += med->pos[l];
		}
		else {
			if (med->pos[l] > L2_max) {
				L2_max = med->pos[l];
			}	
			if (med->pos[l] < L2_min) {
				L2_min = med->pos[l];
			}
			C2 += med->pos[l];
		}
	}
	C1 /= med->n;  //Normalize the center of mass.
	C2 /= med->ntot - med->n;
	printf("[TEST] Boundaries: L1_min = %f, L1_max = %f, and L2_min = %f, L2_max = %f.\n", L1_min, L1_max, L2_min, L2_max);
	printf("[TEST] Barycenter: C1 = %f, and C2 = %f. L1_max/C1 = %f, and L2_max/C2 = %f.\n", C1, C2, L1_max/C1, L2_max/C2);
}

/**
 * Prints the positions of the atoms in the given medium to an ASCII xyz file of given filename. Each line contains the coordinates of an atom.
 * @note Use for checking purpose only.
 */
void save_medium_pos(Medium* med, const char* fname) {
	printf("[INFO] Printing the atomic positions to file '%s'...\n", fname);
	FILE* outfp = fopen(fname, "w");
	int i, c; //Atom index, coordinate.
	fprintf(outfp, "%d\t%d\n", med->d, med->n);
	for (i = 0; i < med->n; i++) { //Loop on atoms.
		for (c = 0; c < med->d; c++) { //Loop on coordinates.
			fprintf(outfp, "\t%.9f", med->pos[i*med->d + c]);
			//if (c != med->d - 1) {
			//	fprintf(outfp, "\t");
			//}
		}
		fprintf(outfp, "\n");
	}
	fclose(outfp);
}

/**
 * Test the computation of the inverse scattering amplitude and the cross section.
 */
void test_medium_invf(Medium* med) {
	print_param_medium(med);
	dcomplex k =  6.;
	print_cross_section(med, k);
	k = 10.;
	print_cross_section(med, k);
	k = 6. - 3.*I;
	print_cross_section(med, k);
}

int main(int argc, char** argv) {
	
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	
	//2D case:
	char* args1[] = {"dimension=2", "natom=100", "model=hardsphere", "alpha=0.15", "shape=ball", "ratio=1.0"};
	parse_medium(med, 6, args1);
	test_medium_invf(med);
	
	//3D case:
	char* args2[] = {"dimension=3", "natom=100", "model=hardsphere", "alpha=0.26", "shape=ball", "ratio=1.0"};
	parse_medium(med, 6, args2);
	test_medium_invf(med);
	
	del_medium(med);
	free(med);
	return 0;
}
