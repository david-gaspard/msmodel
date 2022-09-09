/****
 * @date Created on 2021-04-01 at 16:16:27 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing some test functions to test the eigenstate utilities.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include "eigenstate_util.h"
#include "medium_util.h"

void test_parse_eigenstate() {
	char *args[] = {
		"seed=17",
		"maxit=20",
		"toler=1e-9",
		"kinit=1-1i 2-3i",
		//"realtime=1.2458e3",
		"conv0=1",
		"path0=1-1i 1e-2-2e-3i 4.2+2.3i 4.1+2.2i 4.05+2.1i 4.02+2.01i 4.0+2.0i",
		"eigv0=0.5+0.3i 1e-2-2e-3i 4+2i",
		"conv1=0",
		"path1=2-3i 1.8-3.4i 1.5-3.7i 1.2-3.8i 1.1-3.9i 1.05-3.95i 1.01-4.002i 4.0-4.0i",
		"eigv1=  0.7+0.2i  -4.1-0.3i  -1.5+2.4i  ",
	};
	int narg = sizeof(args)/sizeof(args[0]);
	
	Eigenstate* eigst = (Eigenstate*)calloc(1, sizeof(Eigenstate));
	parse_eigenstate(eigst, narg, args);
	
	if (eigst->data[0].eigv[0] == 0.5+0.3*I && eigst->data[1].eigv[0] == 0.7+0.2*I) {
		printf("[INFO] OK, the eigv[0] values are correct...\n");
	}
	if (eigst->data[0].len == 3 && eigst->data[1].len == 3) {
		printf("[INFO] OK, the 'len' values are correct...\n");
	}
	if (eigst->data[0].path[0] == 1.-1.*I && eigst->data[1].path[0] == 2.-3.*I) {
		printf("[INFO] OK, the path[0] values are correct...\n");
	}
	if (eigst->data[0].niter == 7 && eigst->data[1].niter == 8) {
		printf("[INFO] OK, the 'niter' values are correct...\n");
	}
	if (eigst->data[0].conv == 1 && eigst->data[1].conv == 0) {
		printf("[INFO] OK, the 'conv' values are correct...\n");
	}
	
	save_eigenstate(eigst, "test_eigenstate.dat");
	del_eigenstate(eigst);
	free(eigst);
}

int main(int argc, char** argv) {
	test_parse_eigenstate();
	return 0;
}
