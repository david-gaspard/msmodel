/****
 * @date Created on 2021-02-26 at 11:09:04 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to test the solution method of the classical diffusion equation.
 ***/
#include <stdio.h>
#include "classical_diffusion_util.h"

/**
 * Tests the computation in Taylor series of the regular solution of the Helmholtz equation.
 */
void test_helmholtz() {
	int d = 3;
	double k = 1.;
	double r = 3.;
	double P = helmholtz_series(d, k, r);
	double drP = dr_helmholtz(d, k, r);
	double dkP = dk_helmholtz(d, k, r);
	double dkdrP = dkdr_helmholtz(d, k, r);
	printf("[INFO] At d=%d, k=%g, r=%g :\n", d, k, r);
	printf("[INFO] P(d,k,r)         = %.16g.\n", P);
	printf("[INFO] D_r P(d,k,r)     = %.16g.\n", drP);
	printf("[INFO] D_k P(d,k,r)     = %.16g.\n", dkP);
	printf("[INFO] D_k D_r P(d,k,r) = %.16g.\n", dkdrP);
}

void test_diffusion_solve() {
	int d = 3;
	double radius = 7.81593;
	double nsigma = 0.349025;
	double beta, kappa, betahyp, kappahyp;
	int conv, convhyp;
	//diffusion_solve(d, radius, nsigma, &beta, &kappa, &betahyp, &kappahyp);
	diffusion_solve_v1(d, radius, nsigma, &beta, &kappa, &conv, &betahyp, &kappahyp, &convhyp);
	printf("[INFO] Std sol:\tbeta  = %.16g\tkappa  = %.16g\t(conv  = %d)\n", beta, kappa, conv);
	printf("[INFO] Hyp sol:\tbetah = %.16g\tkappah = %.16g\t(convh = %d)\n", betahyp, kappahyp, convhyp);
}

int main(int argc, char** argv) {
	//test_helmholtz();
	test_diffusion_solve();
	return 0;
}
