/****
 * @date Created on 2020-08-13 at 15:47:28 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to test the complex map object, creation and manipulation tools.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include "medium_util.h"
#include "complex_map_util.h"
#include "config_io.h"

int main(int argc, char** argv) {
	
	double xmin =  0.0;
	double xmax =  4.0;
	double ymin = -3.0;
	double ymax =  0.0;
	int nx = 5;
	int ny = 4;
	int nseed = 1;
	
	printf("[INFO] Creating complex map 'cmap'... "); fflush(stdout);
	ComplexMap* cmap = (ComplexMap*)calloc(1, sizeof(ComplexMap));
	init_complexmap(cmap, nx, ny, xmin, xmax, ymin, ymax, nseed, "Testing purpose only");
	printf("Done\n"); fflush(stdout);
	
	printf("[TEST] Testing the argument of the complex map:\n");
	dcomplex z;
	int l;
	for (l = 0; l < cmap->ntot; l++) {
		z = cmap->get_arg(cmap, l);
		printf("\t%f%+fi", creal(z), cimag(z));
		if ((l+1)%nx == 0) {
			printf("\n");
		}
	}
	
	printf("[TEST] Sample map:\n");
	set_value_cmap(cmap, 0, 0,   1.39415823852+5.1242*I);
	set_value_cmap(cmap, 1, 0,   2.39415823852+1.0254*I);
	set_value_cmap(cmap, 2, 0,   3.39415823852-2.5145*I);
	set_value_cmap(cmap, 3, 0,   4.39415823852-0.0124*I);
	set_value_cmap(cmap, 4, 0,   5.39415823852+2.9514*I);
	
	set_value_cmap(cmap, 0, 1,   0.25415823346-0.1249*I);
	set_value_cmap(cmap, 1, 1,   2.25415823346-1.0259*I);
	set_value_cmap(cmap, 2, 1,  10.25415823346+3.5149*I);
	set_value_cmap(cmap, 3, 1,  30.25415823346+4.0129*I);
	set_value_cmap(cmap, 4, 1,  68.25415823346-1.9519*I);
	
	set_value_cmap(cmap, 0, 2,  -3.192453258725+2.1242*I);
	set_value_cmap(cmap, 1, 2,   0.192453258725-1.0254*I);
	set_value_cmap(cmap, 2, 2,  15.192453258725+0.5145*I);
	set_value_cmap(cmap, 3, 2,  54.192453258725-1.0124*I);
	set_value_cmap(cmap, 4, 2, 129.192453258725+4.9514*I);
	
	set_value_cmap(cmap, 0, 3,  -8.92453258435-5.1243*I);
	set_value_cmap(cmap, 1, 3,  -4.92453258435-1.0253*I);
	set_value_cmap(cmap, 2, 3,  18.92453258435+2.5143*I);
	set_value_cmap(cmap, 3, 3,  76.92453258435-0.0123*I);
	set_value_cmap(cmap, 4, 3, 188.92453258435+2.9513*I);
	
	print_complexmap(cmap, 5);
	
	double minreal, maxreal;
	find_min_max_cmap(cmap, 1);
	
	printf("[TEST] Mininum real = %f. Maximum real = %f\n", cmap->colrule->hmin, cmap->colrule->hmax);
	
	//Test bilinear interpolation:
	z = -0.2 + 0.5*I;
	dcomplex fz = interpolate_bilinear(cmap, z);
	printf("[TEST] Interpolated value at z = %f%+fi is f(z) = %f%+fi.\n", creal(z), cimag(z), creal(fz), cimag(fz));
	
	//char* fname = "cmap_test"; //Without extension.
	//plot_complexmap_real(cmap, fname);
	
	//printf("[TEST] Laplace map:\n");
	//ComplexMap* laplmap = new_laplacemap(cmap, "Testing Laplace map");
	//print_complexmap(laplmap, 5);
	//fname = "laplmap_test"; //Without extension.
	//plot_complexmap_real(laplmap, fname);
	
	//char* fname = "cmap_test.dat";
	//printf("[TEST] Complex map and metadata output to file '%s'...\n", fname);
	//int d = 2;
	//int n = 150;
	//double ratio = 2.5;
	//double alpha = 1.3;
	//uint64_t seed = 14;
	//cmap->realtime = 1659.23;
	//Medium* med = new_medium(d, n, ratio, alpha);
	//fill_medium(med, seed);
	//save_data(cmap, med, fname);
	
	const char* fname_in = "sample.conf";
	Config* conf = new_config(fname_in);
	
	printf("[INPUT] Dimension = %d\n", conf->med->d);
	printf("[INPUT] Nb. Atoms = %d\n", conf->med->n);
	printf("[INPUT] Nb. Seeds = %d\n", conf->cmap->nseed);
	printf("[INPUT] Map Title = '%s'\n", conf->cmap->title);
	printf("[INPUT] Real Time = %g s\n", conf->cmap->realtime);
	printf("[INPUT] Rek Range = %g:%g\n", conf->cmap->xmin, conf->cmap->xmax);
	printf("[INPUT] Imk Range = %g:%g\n", conf->cmap->ymin, conf->cmap->ymax);
	printf("[INPUT] RekSample = %d\n", conf->cmap->nx);
	printf("[INPUT] ImkSample = %d\n", conf->cmap->ny);
	printf("[INPUT] CMap Data = \n");
	print_complexmap(conf->cmap, 5);
	printf("[TEST] Compare with original data = \n");
	print_complexmap(cmap, 5);
	
	del_config(conf);
	del_complexmap(cmap);
	free(cmap);
	
	return 0;
}
