/****
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created on 2020-07-26 at 12:46 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Main file of the Multiple Scattering Program MSModel.
 ***/
#include <stdlib.h>         /* Standard Library providing exit() function */
#include <stdio.h>          /* Standard Library for Input and Output */
#include <unistd.h>         /* Library providing the UNIX Getopt Utility */
#include "common_util.h"    /* Import some General Purpose Functions such as the UUID Generator */
#include "ms_core.h"        /* Import the core functions of the multiple scattering program */
/**
 * Number of OpenMP threads which should be equal to the number of cores, rather than threads, to avoid huge load imbalancing
 * after each random geometry due to hyper-threading saturation (see: $ cat /proc/cpuinfo | grep -wi "cores\|siblings").
 * Loss of performance can be expected beyond the number of cores.
 */
static const int  NUMTHREADS = 4; //Default number of OpenMP threads.
static const char OUTPUTDIR[] = "out/";  //Name of the output directory with trailing slash.
static const char GREETING[] = "====== This is MSModel ======\n";
static const char USAGE[] = "[USAGE] %s [-p nthread] file1.conf [file2.conf file3.conf ...]\n";

/**
 * Main function of the Multiple Scattering Program MSModel.
 */
int main(int argc, char **argv) {
	
	printf(GREETING);
	int nthread = NUMTHREADS;  //Default thread number used by OpenMP.
	int opt = 0;
	while ((opt = getopt(argc, argv, "?p:")) != -1) {//Loop on option flags.
		switch (opt) {//NB: Other errors are treated in the corresponding object constructors.
			case 'p':
				if (sscanf(optarg, "%d", &nthread) != 1 || nthread <= 0) {
					printf("[ERROR] Invalid thread number '%s', aborting...\n", optarg);
					exit(1);
				}
				break;
			case '?':
				printf(USAGE, argv[0]);
				exit(0);
				break;
		}
	}
	if (optind == argc) {//If no further arguments.
		printf("[ERROR] No input file...\n");
		printf(USAGE, argv[0]);
		exit(1);
	}
	char* infname;
	char outfname[120];  //Output data file name.
	int i;
	for (i = optind; i < argc; i++) {//Loop on input files.
		infname = argv[i];
		if (no_file(infname)) {
			printf("[ERROR] No file '%s' found, aborting...\n", infname);
			exit(1);
		}
		if (no_text_file(infname)) {
			printf("[ERROR] File '%s' is not a text/ascii file, aborting...\n", infname);
			exit(1);
		}
		Config* conf = new_config(infname);    //Load and parse all the data.
		generate_filename(conf, OUTPUTDIR, outfname); //Generate the filename according to the input file configuration "conf".
		ensure_path(outfname);                 //Ensure or create a local directory for the output file.
		prepare_file(outfname, "#", "Computed");   //Prepare the output data file with a header timestamp. Removes any existing file of the same name (very unlikely).
		ms_core_main(conf, outfname, nthread); //Performs all the computations planned in the given file.
		del_config(conf);
	}
	return 0;  //or exit(EXIT_SUCCESS);
}
