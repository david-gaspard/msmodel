/****
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created on 2020-07-25 at 12:33
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Main file to test the output of the thread safe pseudo-random generator.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include "nru_random.h"   //Thread safe generator.

/**
 * Main function returning a list of random numbers.
 */
int main(int argc, char** argv) {
	
	uint64_t nseed = 5;  //Number of seeds.
	size_t nval = 20000;    //Number of computed random values for a given seed.
	
	double* table = (double*)calloc(nseed*nval, sizeof(double));  //OMP-Shared list of random values.
	
	#pragma omp parallel shared(table)
	{
		Urandom rnd;  //OMP-Private random generator (thread safe).
		uint64_t s;
		size_t i;
		#pragma omp for
		for (s = 1; s <= nseed; s++) {//Parallelized loop on seeds.
			init_random(&rnd, s); //Initialize with the given seed.
			for (i = 0; i < nval; i++) {
				//table[(s-1)*nval+i] = random_double(&rnd);
				table[(s-1)*nval+i] = random_normal(&rnd);
			}
		}
	}
	//char filename[50];
	//sprintf(filename, "test_randomseq_s%lu.tab", seed);
	//FILE* outfp = fopen(filename, "w");
	
	FILE* outfp = fopen("test_random_normal.tab", "w");
	uint64_t s, i;
	for (s = 1; s <= nseed; s++) {
		fprintf(outfp, "SEED=%lu\t\t", s);
	}
	fprintf(outfp, "\n");
	for (i = 0; i < nval; i++) {
		for (s = 1; s <= nseed; s++) {
			fprintf(outfp, "%.12f\t", table[(s-1)*nval+i]);
		}
		fprintf(outfp, "\n");
	}
	fclose(outfp);
	
	free(table);
	return 0;
}
