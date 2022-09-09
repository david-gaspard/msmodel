/****
 * @date Created on 2020-09-07 at 17:07:28 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code to test the configuration parser.
 ***/
#include <stdio.h>
#include <string.h>
#include "common_util.h"
#include "medium_util.h"
#include "complex_map_util.h"
#include "chistogram_util.h"
#include "config_io.h"

int main(int argc, char** argv) {
	if (argc > 1) {
		Config* conf = new_config(argv[1]);
		print_param_medium(conf->med);
		char stripfname[strlen(argv[1])], fname[strlen(argv[1])+8];
		replace_ext(argv[1], "", stripfname);
		sprintf(fname, "test_%s", stripfname);
		int i1, i2;
		for (i1 = 0; i1 < conf->ncmap; i1++) {
			print_param_complexmap(conf->cmap+i1, 1);
			print_complexmap(conf->cmap, 5);
			export_tikz_map(conf->cmap+i1, i1, fname);
		}
		for (i2 = 0; i2 < conf->nchist; i2++) {
			print_param_chistogram(conf->chist+i2, 1);
			export_tikz_chistogram(conf->chist+i2, i1+i2, fname);
		}
		printf("[INFO] Now trying to free data...\n"); fflush(stdout);
		del_config(conf);
	}
	else {
		printf("[USAGE] %s file.cfg \n", argv[0]);
	}
	return 0;
}
