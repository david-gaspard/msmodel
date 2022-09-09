/****
 * @date Created on 2021-03-12 at 10:51:49 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the testing function for the resonance root-finding in the square well
 * effective approximation model of the multiple scattering in a ball-shaped random Lorentz gas.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "green_bessel.h"
#include "domain_util.h"
#include "medium_util.h"
#include "square_well_resonances.h"

void test_pade_roots() {
	dcomplex pade[] = {13., -2., 4., 9., 3., -5., 6.};
	int i, nz = sizeof(pade)/sizeof(dcomplex) - 1;
	printf("[INFO] Polynomial: P(x) = ");
	for (i = 0; i <= nz; i++) {
		printf("(%g%+gi)x^%d", creal(pade[i]), cimag(pade[i]), i);
		if (i != nz)
			printf(" + ");
	}
	printf("\n");
	dcomplex roots[nz];
	pade_roots(nz, pade, roots);
	for (i = 0; i < nz; i++) {
		printf("[INFO] Root #%d = %.12g %+.12gi \n", i, creal(roots[i]), cimag(roots[i]));
	}
}

void print_carray(int n, dcomplex* z, char* msg) {
	printf("[INFO] %s zs = {", msg);
	int i;
	for (i = 0; i < n; i++) {
		printf("%g%+gi", creal(z[i]), cimag(z[i]));
		if (i != n-1)
			printf(", ");
	}
	printf("}, size=%d.\n", n);
}

void test_exclude_from_domain() {
	dcomplex zs[] = {1+I, 13.+8*I, -2.-I, 4.+3*I, 9.+10*I, 3.-7*I, -5.-4*I, 6.+2.*I};
	int nz = sizeof(zs)/sizeof(dcomplex);
	Domain dom = {//Rectangle[{-3, -6}, {10, 7}]
		.xmin = -3.,
		.xmax = 10.,
		.ymin = -8.,
		.ymax =  7.
	};
	print_carray(nz, zs, "Init");
	exclude_from_domain(&nz, zs, &dom);
	print_carray(nz, zs, "Excluded");
}

void test_sort_carray() {
	dcomplex zs[] = {13.+8*I, -2.-I, 4.+3*I, 9.+10*I, 3.-7*I, -5.-4*I, 6.+2.*I};
	int nz = sizeof(zs)/sizeof(dcomplex);
	sort_carray(nz, zs, creal);
	print_carray(nz, zs, "Real Sorted");
	sort_carray(nz, zs, cabs);
	print_carray(nz, zs, "Abs  Sorted");
}

void test_uniq_carray() {
	dcomplex zs[] = {13.+8*I, 13.+8*I, -2.-I, -2-I, -2-I, 4.+3*I, 9.+10*I, 9+10*I, 3.-7*I, -5.-4*I, -5-4*I, -5-4*I, 6.+2.*I};
	int nz = sizeof(zs)/sizeof(dcomplex);
	print_carray(nz, zs, "Init");
	uniq_carray(&nz, zs, 1e-9);
	print_carray(nz, zs, "Uniq");
}

//void bessel_k_ratio_function(ResonanceEquation* req, dcomplex k, dcomplex* f) {
//	int twonu = 2*req->l + req->med->d - 2;  //Two times nu = l + (d-2)/2.
//	*f = bessel_k_ratio(twonu, -I*k*radius_of_ball(req->med));
//}

void test_solve_secant() {
	
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {"dimension=3", "natom=100", "model=hardsphere", "alpha=0.1", "shape=ball", "ratio=1"};
	parse_medium(med, 6, args);
	Domain dom = {
		.xmin = -0.5,
		.xmax =  5.0,
		.ymin = -3.5,
		.ymax =  0.5
	};
	ResonanceEquation* req = (ResonanceEquation*)calloc(1, sizeof(ResonanceEquation));
	req->med = med;
	req->dom = &dom;
	req->function = resonance_function;  //bessel_k_ratio_function;
	req->l = 5;
	dcomplex root = 0.7 - 1.5*I;  //Initial guess.
	if (solve_secant(req, &root)) //Run the secant iteration.
		printf("[INFO] Successfully converged to root = %.15g %+.15gi sp^-1.\n", creal(root), cimag(root));
	else
		printf("[WARN] Not converged. Result = %.15g %+.15gi sp^-1.\n", creal(root), cimag(root));
	
	del_medium(med);
	free(med);
	free(req);
}

void print_pade(int nz, dcomplex* pade) {
	int i;
	printf("[INFO] Pade P(x) = [");
	for (i = 0; i <= nz; i++) {
		printf("(%g%+gi)x^%d", creal(pade[i]), cimag(pade[i]), i);
		if (i != nz)
			printf(" + ");
	}
	printf("]/[");
	for (i = 0; i <= nz; i++) {
		if (i != 0)
			printf("(%g%+gi)x^%d", creal(pade[i+nz+1]), cimag(pade[i+nz+1]), i);
		else
			printf("1");
		if (i != nz)
			printf(" + ");
	}
	printf("]\n");
}

void test_find_resonances() {
	
	Medium* med = (Medium*)calloc(1, sizeof(Medium));
	char* args[] = {"dimension=2", "natom=100", "model=hardsphere", "alpha=0.1", "shape=ball", "ratio=1"};
	parse_medium(med, 6, args);
	Domain dom = {//Domain of root-finding.
		.xmin = -0.5,
		.xmax =  5.0,
		.ymin = -3.5,
		.ymax =  0.5
	};
	ResonanceEquation* req = (ResonanceEquation*)calloc(1, sizeof(ResonanceEquation));
	req->med = med;
	req->dom = &dom;
	req->function = resonance_function;  //bessel_k_ratio_function;
	req->l = 5;
	
	int nz = (int)ceil(length_x(&dom)*radius_of_ball(med) + 1);  //Roughly estimated number of resonances with margin for all partial waves.	
	int i, np = 2*nz + 1;  //required number of points.
	dcomplex* pade = (dcomplex*)calloc(np, sizeof(dcomplex));  //Allocate space for the Pade approximant, stored as [a0, a1, a2, ..., an, b1, b2, ..., bn], with b0=1 by convention.
	dcomplex* zs   = (dcomplex*)calloc(np, sizeof(dcomplex));  //Allocate space for the arguments of the sample points in the complex domain.
	dcomplex* fs   = (dcomplex*)calloc(np, sizeof(dcomplex));  //Allocate space for the values of the function at the sample points.
	for (i = 0; i < np; i++) {//#2: Take random uniform samples of the resonance equation in the k-plane.
		zs[i] = halton_point_2_3(req->dom, i+1); //Pick sample point using low-discrepancy Halton-van der Corput's (2,3) sequence.
		req->function(req, zs[i], &fs[i]);       //Compute the function at the given point k = zs[i].
		//printf("[INFO] Sample point z[%d] = %g%+gi \n", i, creal(zs[i]), cimag(zs[i]));
	}
	pade_fit(nz, pade, zs, fs);  //#3: Fit the samples with Pade approximants by solving the linear system (calling LAPACK).
	print_pade(nz, pade);
	
	dcomplex* roots = (dcomplex*)calloc(nz, sizeof(dcomplex));
	pade_roots(nz, pade, roots); //#4: Root-find the numerator of the Pade approximant with the companion matrix method (with LAPACK).
	for (i = 0; i < nz; i++) {
		printf("[INFO] Root #%d = %.15g%+.15gi \n", i+1, creal(roots[i]), cimag(roots[i]));
	}
	printf("[INFO] Filtering roots...\n");
	exclude_from_domain(&nz, roots, &dom); //#7: Again exclude the roots outside of the domain "dom", since secant steps could have moved the roots out of the domain.
	for (i = 0; i < nz; i++) {
		printf("[INFO] Root #%d = %.15g%+.15gi \n", i+1, creal(roots[i]), cimag(roots[i]));
	}
	printf("[INFO] Polishing roots...\n");
	for (i = nz-1; i >= 0; i--) {//#6: Polish the root estimates with the secant method, dealing with one root at a time.
		if (!solve_secant(req, &roots[i])) {//If secant method not converged, then put the root at the end of the array, and reduce the logical array size.
			nz--;  //Reduce logical size of the array.
			cswap(&roots[i], &roots[nz]);
		}
	}
	exclude_from_domain(&nz, roots, &dom); //#7: Again exclude the roots outside of the domain "dom", since secant steps could have moved the roots out of the domain.
	sort_carray(nz, roots, creal); //#8: Sort the roots in ascending order of the real part.
	uniq_carray(&nz, roots, 1e-9); //#9: Remove duplicate adjacent roots within small tolerance, still larger than the machine epsilon.
	for (i = 0; i < nz; i++) {
		printf("[INFO] Root #%d = %.15g%+.15gi \n", i+1, creal(roots[i]), cimag(roots[i]));
	}
	free(roots);
	free(pade);
	free(zs);
	free(fs);
}

int main(int argc, char** argv) {
	
	//test_pade_roots();
	//test_exclude_from_domain();
	//test_sort_carray();
	//test_uniq_carray();
	//test_solve_secant();
	test_pade();
	
	return 0;
}
