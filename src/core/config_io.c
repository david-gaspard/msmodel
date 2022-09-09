/****
 * @date Created on 2020-09-07 at 13:00:12 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing the parser functions to import the settings and data from a configuration file *.conf or *.cfg.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "common_util.h"
#include "config_type.h"
#include "medium_util.h"
#include "complex_map_util.h"
#include "complex_cut_util.h"
#include "log_cut_util.h"
#include "chistogram_util.h"
#include "eigenstate_list_util.h"
#include "random_function_util.h"
#include "diff_cross_section_util.h"
#include "imag_mu_histogram_util.h"
#include "wave_function_util.h"
#include "wave_function_cut_util.h"

/**
 * Strip the trailing white spaces at the end of the given "line".
 * Also strip the spaces surrounding the key/value separator "=" or ":", if any.
 */
void strip_line(char* line) {
	int i, x = 0;     //x=target index of the character to be written.
	int valflag = 0;  //Flag becomes 1 after equal sign, and 2 within value field.
	for (i = 0; line[i]; i++) {
		if (!is_blank(line[i]) || valflag == 2) {
			line[x] = line[i];
			x++;
			if (line[i] == '=' || line[i] == ':') {//Encounters the key/value separator.
				valflag = 1;
			}
			else if (valflag == 1) {//Really enters the field zone.
				valflag = 2;
			}
		}
	}
	line[x] = '\0';
	i = x-1;
	while (i >= 0 && is_blank(line[i])) {//Strip trailing white spaces.
		line[i] = '\0';
		i--;
	}
}

/**
 * Initializes the given "lines" while running through the file "fname".
 * The "mode" can be either allocation of whole "lines" (mode=0), allocation of individual lines (mode=1), or reading (mode=2).
 * Returns the number of lines found in file "fname".
 */
int init_lines(const char* fname, char*** lines, int mode) {
	FILE* fp = fopen(fname, "r");
	int l = 0;        //Current line index.
	if (fp == NULL) { //Failed to open file.
		return l;
	}
	int i = 0;        //Index of current character in a line. Also number of characters in the current line (excluding comments).
	int comment = 0;  //Flag is 1 in comment.
	char c;           //Current character.
	while ((c = fgetc(fp)) != EOF) {
		if (is_blank(c) && i == 0) {//Skip space at beginning of line.
			continue;
		}
		else if (c == '#' || c == '!' || c == ';') {//Activate comment flag.
			comment = 1;
		}
		else if (c == '\n') {
			if (i > 0) {
				if (mode == 1) {//If individual line allocation mode.
					//printf("Line[%2d] allocating %d characters.\n", l, i+1);
					(*lines)[l] = (char*)calloc(i+1, sizeof(char));
				}
				else if (mode == 2) {//Appends null character.
					(*lines)[l][i] = '\0';
					strip_line((*lines)[l]); //Strip the line.
				}
				l++;     //Go to next line.
				i = 0;   //Reset current number of in-line characters to zero.
			}
			comment = 0; //We are back out of comment.
		}
		else if (!comment) {//Read only uncommented characters.
			if (mode == 2) {//If read mode.
				(*lines)[l][i] = c;
			}
			i++;
		}
	}
	if (mode == 0) {//If lines allocation mode.
		//printf("[INFO] Allocating %d lines.\n", l);
		*lines = (char**)calloc(l, sizeof(char*));
	}
	fclose(fp);
	return l;
}

/**
 * Reads the lines in file "fname" and store them in "*lines". The count of lines is stored in "*lc".
 * Omit the comments and the blank space at the beginning of each line.
 */
void read_lines(const char* fname, int* lc, char*** lines) {
	*lc = init_lines(fname, lines, 0);  //Allocate memory for the whole "lines", counting the number of lines in "fname".
	*lc = init_lines(fname, lines, 1);  //Allocate memory for each individual line, coutning the number of characters in each line.
	*lc = init_lines(fname, lines, 2);  //Read and store each line. Stripping useless white spaces.
}

/**
 * Removes the given lines from memory. Do not forget to call after "init_lines()".
 */
void del_lines(int lc, char** lines) {
	int l;
	for (l = 0; l < lc; l++) {
		free(lines[l]);
	}
	free(lines);
}

/**
 * Prints lines to stdout for testing purposes.
 */
void print_lines(int lc, char** lines) {
	int l;
	for (l = 0; l < lc; l++) {
		printf("Line[%2d] = '%s'\n", l, lines[l]);
	}
}

/**
 * Prints lines to stdout for testing purposes (one-liner version).
 */
void print_lines_short(int lc, char** lines) {
	int l;
	printf("{");
	for (l = 0; l < lc; l++) {
		printf("'%s'", lines[l]);
		if (l != lc-1) {
			printf(", ");
		}
	}
	printf("}\n");
}

/**
 * Count the number of arguments in "lines" before the first config section encountered and starting from index "l0".
 * The total line count in "lines" is "lc".
 */
int count_args(int l0, int lc, char** lines) {
	int l;
	for (l = l0; l < lc; l++) {
		if (lines[l][0] == '[') {//Encountered a config section.
			return l-l0;
		}
	}
	return lc-l0; //If no config section encountered.
}

/**
 * Initializes the context, i.e., the many items involved in the numerical simulations, from the given instruction "lines".
 * The "mode" can be either memory allocation (mode=0), or parsing (mode=1).
 */
void init_config(Config* conf, int lc, char** lines, int mode) {
	conf->nmed = 0;   //Number of medium, should be 1.
	conf->ncmap = 0;  //Number of k-plane complex map.
	conf->nccut = 0;  //Number of k-plane complex cut.
	conf->nlcut = 0;  //Number of k-plane log cut.
	conf->nchist = 0; //Number of mu-plane eigenvalue histogram.
	conf->neigls = 0; //Number of eigenstates list objects.
	conf->nrfun = 0;  //Number of random function objects.
	conf->ndcs = 0;   //Number of differential cross section objects.
	conf->nihist = 0; //Number of Im(mu) histogram objects.
	conf->nwfun = 0;  //Number of wave function objects.
	conf->nwfcut = 0; //Number of wave function cut objects.
	int l, narg;
	for (l = 0; l < lc; l++) {//Loop to count the number of planned simulation items.
		if (strcmp(lines[l], "[medium]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					//printf("[INFO] Type = '%s', Arg count = %d\n", lines[l], narg);
					parse_medium(conf->med, narg, lines+l+1); //Parse the medium with the given arguments.
				}
				conf->nmed++;
			}
		}
		else if (strcmp(lines[l], "[kplane]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_complexmap(conf->cmap+conf->ncmap, narg, lines+l+1); //Parse the complex map with the given arguments.
				}
				conf->ncmap++;
			}
		}
		else if (strcmp(lines[l], "[kcut]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_complexcut(conf->ccut+conf->nccut, narg, lines+l+1); //Parse the complex cut with the given arguments.
				}
				conf->nccut++;
			}
		}
		else if (strcmp(lines[l], "[logkcut]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_logcut(conf->lcut+conf->nlcut, narg, lines+l+1); //Parse the fast cut with the given arguments.
				}
				conf->nlcut++;
			}
		}
		else if (strcmp(lines[l], "[muplane]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					//printf("[INFO] Type = '%s', Arg count = %d\n", lines[l], narg);
					parse_chistogram(conf->chist+conf->nchist, narg, lines+l+1); //Parse the complex histogram with the given arguments.
				}
				conf->nchist++;
			}
		}
		else if (strcmp(lines[l], "[eigenstate_list]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					if (conf->nmed == 0) {//Check for existence of medium.
						printf("[ERROR] Cannot create eigenstates without a medium, aborting...\n");
						exit(EXIT_FAILURE);
					}
					//printf("[INFO] Type = '%s', Arg count = %d\n", lines[l], narg);
					parse_eigenstate_list(conf->eigls+conf->neigls, conf->med, narg, lines+l+1); //Parse the eigenstates list object with the given arguments.
				}
				conf->neigls++;
			}
		}
		else if (strcmp(lines[l], "[random_function]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_random_function(conf->rfun+conf->nrfun, narg, lines+l+1); //Parse the random function with the given arguments.
				}
				conf->nrfun++;
			}
		}
		else if (strcmp(lines[l], "[diff_cross_section]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_diff_cross_section(conf->dcs+conf->ndcs, narg, lines+l+1); //Parse the random function with the given arguments.
				}
				conf->ndcs++;
			}
		}
		else if (strcmp(lines[l], "[imag_mu]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_ihist(conf->ihist+conf->nihist, narg, lines+l+1); //Parse the random function with the given arguments.
				}
				conf->nihist++;
			}
		}
		else if (strcmp(lines[l], "[wavefunction]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_wfun(conf->wfun+conf->nwfun, narg, lines+l+1); //Parse the random function with the given arguments.
				}
				conf->nwfun++;
			}
		}
		else if (strcmp(lines[l], "[wavefunction_cut]") == 0) {
			narg = count_args(l+1, lc, lines);
			if (narg > 0) {//Check for arguments.
				if (mode == 1) {
					parse_wfcut(conf->wfcut+conf->nwfcut, narg, lines+l+1); //Parse the random function with the given arguments.
				}
				conf->nwfcut++;
			}
		}
	}
	if (mode == 0) {//Allocate space for simulation items.
		if (conf->nmed != 1) {
			if (conf->nmed == 0) {
				printf("[ERROR] No medium found, aborting...\n");
			}
			else {
				printf("[ERROR] Found too many media (%d), aborting...\n", conf->nmed);
			}
			exit(EXIT_FAILURE);
		}
		//printf("[INFO] Allocating one medium...\n");
		conf->med = (Medium*)calloc(1, sizeof(Medium));
		//printf("[INFO] Allocating %d k-plane...\n", conf->ncmap);
		conf->cmap = (ComplexMap*)calloc(conf->ncmap, sizeof(ComplexMap));
		//printf("[INFO] Allocating %d k-cut...\n", conf->nccut);
		conf->ccut = (ComplexCut*)calloc(conf->nccut, sizeof(ComplexCut));
		//printf("[INFO] Allocating %d log k-cut...\n", conf->nlcut);
		conf->lcut = (LogCut*)calloc(conf->nlcut, sizeof(LogCut));
		//printf("[INFO] Allocating %d mu-plane...\n", conf->nchist);
		conf->chist = (Chistogram*)calloc(conf->nchist, sizeof(Chistogram));
		//printf("[INFO] Allocating %d eigenstates list...\n", conf->neigls);
		conf->eigls = (EigenstateList*)calloc(conf->neigls, sizeof(EigenstateList));
		//printf("[INFO] Allocating %d random functions...\n", conf->nrfun);
		conf->rfun = (RandomFunction*)calloc(conf->nrfun, sizeof(RandomFunction));
		//printf("[INFO] Allocating %d differential cross sections...\n", conf->ndcs);
		conf->dcs = (DiffCrossSection*)calloc(conf->ndcs, sizeof(DiffCrossSection));
		//printf("[INFO] Allocating %d Im(mu) histograms...\n", conf->nihist);
		conf->ihist = (ImagMuHistogram*)calloc(conf->nihist, sizeof(ImagMuHistogram));
		//printf("[INFO] Allocating %d wave functions...\n", conf->nwfun);
		conf->wfun = (WaveFunction*)calloc(conf->nwfun, sizeof(WaveFunction));
		//printf("[INFO] Allocating %d wave function cuts...\n", conf->nwfcut);
		conf->wfcut = (WaveFunctionCut*)calloc(conf->nwfcut, sizeof(WaveFunctionCut));
	}
}

/**
 * Creates a new configuration structure containing the simulation data.
 * Reads and parses the given configuration file "fname".
 */
Config* new_config(const char* fname) {
	//First step, separate each line of input stripping out the comments.
	//Methodology: First, count number of lines, then allocate "config". Second, runs through each line, count the number of characters, and allocate "config[i]". Third runs through each character and store in config[i].
	int lc;
	char** lines = NULL;
	read_lines(fname, &lc, &lines);
	if (lines == NULL) {
		printf("[ERROR] Failed to read file '%s'...\n", fname);
		exit(EXIT_FAILURE);
	}
	else if (lc == 0) {
		printf("[ERROR] No data found in '%s'...\n", fname);
		exit(EXIT_FAILURE);
	}
	//print_lines(lc, lines);
	//Second step, allocate/initialize simulation items.
	Config* conf = (Config*)calloc(1, sizeof(Config));
	init_config(conf, lc, lines, 0); //Allocate space for simulation items.
	init_config(conf, lc, lines, 1); //Parse simulation items.
	del_lines(lc, lines); //Delete lines and simulation items with del_* functions:
	return conf;
}

/**
 * Deletes the given configuration structure and the items herein.
 * This operation will crash if one of the simulation items (med/cmap/chist/...) is not properly built.
 */
void del_config(Config* conf) {
	del_medium(conf->med); //Deletes the content of the simulation items (med/cmap/chist/...), not the pointers.
	int i;
	for (i = 0; i < conf->ncmap; i++) {
		del_complexmap(conf->cmap+i);
	}
	for (i = 0; i < conf->nccut; i++) {
		del_complexcut(conf->ccut+i);
	}
	for (i = 0; i < conf->nlcut; i++) {
		del_logcut(conf->lcut+i);
	}
	for (i = 0; i < conf->nchist; i++) {
		del_chistogram(conf->chist+i);
	}
	for (i = 0; i < conf->neigls; i++) {
		del_eigenstate_list(conf->eigls+i);
	}
	for (i = 0; i < conf->nrfun; i++) {
		del_random_function(conf->rfun+i);
	}
	for (i = 0; i < conf->ndcs; i++) {
		del_diff_cross_section(conf->dcs+i);
	}
	for (i = 0; i < conf->nihist; i++) {
		del_ihist(conf->ihist+i);
	}
	for (i = 0; i < conf->nwfun; i++) {
		del_wfun(conf->wfun+i);
	}
	for (i = 0; i < conf->nwfcut; i++) {
		del_wfcut(conf->wfcut+i);
	}
	//Then, deletes the pointers themselves:
	//printf("[INFO] Freeing the med/cmap/chist pointers themselves...\n"); fflush(stdout);
	free(conf->med);
	free(conf->cmap);
	free(conf->ccut);
	free(conf->lcut);
	free(conf->chist);
	free(conf->eigls);
	free(conf->rfun);
	free(conf->dcs);
	free(conf->ihist);
	free(conf->wfun);
	free(conf->wfcut);
	free(conf);
}

/**
 * Generates a unique data filename, "fname", for the given configuration "conf" within the given directory "outputdir" (with trailing slash).
 * The string "fname" should be initialized with enough characters (about 100).
 */
void generate_filename(Config* conf, const char* outputdir, char* fname) {
	char uuid[20];         //UUID string, no buffer overflow expected before 10^23 years.
	generate_uuid(uuid);   //Generate a UUID once.
	sprintf(fname, "%sd%d/%s/%s/n%d_%s.dat", outputdir, conf->med->d, shape_to_string(conf->med->shape), conf->med->scmodel->dir, conf->med->n, uuid);
}
