/****
 * @date Created on 2020-08-18 at 17:15:11 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Main file of the graphical and post-treatment sub-program of the multi-scattering model program.
 ***/
#include <stdlib.h>                   /* Standard Library for Exit Commands */
#include <stdio.h>                    /* Standard Library of Input and Output */
#include <string.h>                   /* Standard Library of String Manipulation */
#include <unistd.h>                   /* Library providing the UNIX Getopt Utility */
#include "common_util.h"              /* Impot the General Purpose Functions */
#include "config_io.h"                /* Import the Input Parsing Functions */
#include "medium_util.h"              /* Import the Medium Utilities */
#include "complex_map_util.h"         /* Import the Complex Map Utilities */
#include "complex_cut_util.h"         /* Import the Complex Cut Utilities */
#include "log_cut_util.h"             /* Import the Log Cut Utilities */
#include "chistogram_util.h"          /* Import the Complex Histogram Utilities */
#include "eigenstate_list_util.h"     /* Import the Eigenstates Utilities */
#include "random_function_util.h"     /* Import the Random Function Utilities */
#include "diff_cross_section_util.h"  /* Import the Differential Cross Section Utilities */
#include "imag_mu_histogram_util.h"   /* Import the Utilities for the Histogram of Im(mu) */
#include "wave_function_util.h"       /* Import the Wave Function Utilities */
#include "wave_function_cut_util.h"   /* Import the Wave Function Cut Utilities */
/**
 * Define some constants:
 */
static const char USAGE[] = "[USAGE] %s [-kclxmr] file1.dat [file2.dat file3.dat ...]\n";
static const char COMPILER[] = "lualatex"; //LaTeX compiler. Use "pdflatex" by default, and "lualatex" to exceed TeX's memory limit.
static const char CROPPER[] = "pdfcrop";

/**
 * Compile the tikz picture file "fname.tikz" into a PDF file calling pdflatex.
 * Note that the filename "fname" should not have any extension.
 */
void render_tikzpicture(const char* fname, int crop) {
	char cmd[600+5*strlen(fname)], dir[strlen(fname)+1], bname[strlen(fname)+1];
	sprintf(cmd, "which %s >/dev/null 2>&1", COMPILER);
	if (system(cmd)) {
		printf("[WARN] Command '%s' does not exist. No PDF output.\n", COMPILER);
		return;
	}
	printf("[INFO] Rendering TikZ picture %s.tikz with '%s'...\n", fname, COMPILER);
	split_path(fname, dir, bname); //Separate the basename from the file path from "fname"
	sprintf(cmd, "%s --interaction=nonstopmode --output-directory=%s -jobname=\"%s\" \"\\documentclass{article}\\usepackage[a4paper,margin=20mm]{geometry}\\usepackage{amsmath,amssymb}\\usepackage{pgfplots}\\pgfplotsset{compat=newest}\\pagestyle{empty}\\begin{document}\\input{\\detokenize{%s.tikz}}\\end{document}\" | grep -C 1 -wi --color=auto \"^!\\|^l\\|error\\|undefined\\|unknown\\|missing\\|runaway\\|misplaced\\|multiply\\|exceeded\\|too\\|doesn't\\|ended\\|extra\\|double\\|forget\\|forgotten\\|overfull\\|underfull\"; rm \"%s.aux\" \"%s.log\"", 
		COMPILER, dir, bname, fname, fname, fname);
	system(cmd);
	char pdfname[strlen(fname)+6];
	sprintf(pdfname, "%s.pdf", fname);
	if (no_file(pdfname)) {
		printf("[WARN] No PDF file has been generated.\n");
		return;
	}
	if (crop) {//Attempting to crop the PDF file:
		sprintf(cmd, "which %s >/dev/null 2>&1", CROPPER);
		if (system(cmd)) {
			printf("[WARN] Command '%s' does not exist. The PDF file will not be cropped.\n", CROPPER);
			return;
		}
		sprintf(cmd, "%s --margins 10 \"%s\" \"%s\" >/dev/null", CROPPER, pdfname, pdfname);
		system(cmd);
	}
}

/**
 * Main function of the MSPlot program.
 */
int main(int argc, char** argv) {
	
	int keep  = 0;  //By default, the temporary files *.png *.tikz are deleted.
	int crop  = 0;  //By default, the final PDF file is not cropped.
	int lflag = 0;  //By default, the Laplacian is not taken (cmap).
	int xflag = 0;  //By default, the vertical cut is not plotted (cmap).
	int thlogflag = 0; //By default, the differential cross section is plotted in linear scale log(theta) (dcs). If tflag=1, then plot "theta" in log scale.
	double x = 0.;  //Default value of x (cmap).
	int lmax = -1;  //By default, no partial wave is considered for the computation resonances of the effective square well approx (cmap).
	uint64_t seed = 0;  //By default, no medium outlook is plotted.
	int opt = 0;
	
	while ((opt = getopt(argc, argv, "?kcltx:m:r:")) != -1) {//Loop on option flags.
		switch (opt) {//NB: Other errors are treated in the corresponding object constructors.
			case 'k':
				keep = 1;
				break;
			case 'c':
				crop = 1;
				break;
			case 'l':
				lflag = 1;
				break;
			case 't':
				thlogflag = 1;
				break;
			case 'x':
				xflag = 1;
				if (sscanf(optarg, "%lg", &x) != 1) {
					printf("[WARN] Invalid abscissa after -x option, aborting...\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 'm':
				if (sscanf(optarg, "%lu", &seed) != 1) {
					printf("[WARN] Invalid seed integer after -m option, aborting...\n");
					exit(EXIT_FAILURE);
				}
				break;
			case 'r':
				if (sscanf(optarg, "%d", &lmax) != 1) {
					printf("[WARN] Invalid orbital momentum number after -r option, aborting...\n");
					exit(EXIT_FAILURE);
				}
				break;
			case '?':
				printf(USAGE, argv[0]);
				exit(EXIT_SUCCESS);
				break;
		}
	}
	if (optind == argc) {
		printf("[ERROR] No input file...\n");
		printf(USAGE, argv[0]);
		exit(1);
	}
	char* infname;
	Config* conf;
	int i, j, j1, j2, j3;
	for (i = optind; i < argc; i++) {//Loop on files.
		infname = argv[i];
		if (no_file(infname)) {
			printf("[ERROR] No file '%s' found, aborting...\n", infname);
			exit(1);
		}
		if (no_text_file(infname)) {
			printf("[ERROR] File '%s' is not a text/ascii file, aborting...\n", infname);
			exit(1);
		}
		conf = new_config(infname);
		char stripfname[strlen(infname)], tikzfname[strlen(infname)+7];
		replace_ext(infname, "", stripfname);       //Strip the extension.
		replace_ext(infname, ".tikz", tikzfname);   //Replace the extension.
		prepare_file(tikzfname, "%", "Generated");      //Add short header timestamp. This step removes any pre-existing *.tikz file.
		export_tikz_medium(conf->med, stripfname, seed);  //Show an outlook of the medium for nonzero seed.
		
		for (j1 = 0; j1 < conf->ncmap; j1++) {
			if (lflag) {//Computes the Laplacian.
				ComplexMap* laplmap = (ComplexMap*)malloc(sizeof(ComplexMap));
				init_laplacemap(conf->cmap, laplmap);
				export_tikz_cmap(laplmap, conf->med, lmax, j1, stripfname);
				if (xflag) {
					export_tikz_cut(laplmap, x, stripfname);
				}
				del_complexmap(laplmap);
				free(laplmap);
			}
			else {
				export_tikz_cmap(conf->cmap+j1, conf->med, lmax, j1, stripfname);
				if (xflag) {
					export_tikz_cut(conf->cmap+j1, x, stripfname);
				}
			}
		}
		for (j2 = 0; j2 < conf->nchist; j2++) {
			export_tikz_chistogram(conf->chist+j2, conf->med, j1+j2, stripfname);
		}
		for (j3 = 0; j3 < conf->nwfun; j3++) {
			export_tikz_wfun(conf->wfun+j3, conf->med, j1+j2+j3, stripfname);
		}
		for (j = 0; j < conf->nccut; j++) {
			export_tikz_ccut(conf->ccut+j, conf->med, lflag, stripfname);
		}
		for (j = 0; j < conf->nlcut; j++) {
			export_tikz_lcut(conf->lcut+j, conf->med, stripfname);
		}
		for (j = 0; j < conf->neigls; j++) {
			export_tikz_eigenstate_list(conf->eigls+j, stripfname);
		}
		for (j = 0; j < conf->nrfun; j++) {
			export_tikz_random_function(conf->rfun+j, conf->med, stripfname);
		}
		for (j = 0; j < conf->ndcs; j++) {
			export_tikz_diff_cross_section(conf->dcs+j, conf->med, thlogflag, stripfname);
		}
		for (j = 0; j < conf->nihist; j++) {
			export_tikz_ihist(conf->ihist+j, conf->med, stripfname);
		}
		for (j = 0; j < conf->nwfcut; j++) {
			export_tikz_wfcut(conf->wfcut+j, conf->med, stripfname);
		}
		render_tikzpicture(stripfname, crop); //Compile the TikZ picture into PDF file.
		
		if (!keep) {//If no "keep" flag, remove the temporary *.png *.tikz files (not the *.pdf).
			char cmd[20+strlen(stripfname)+strlen(tikzfname)];
			if (conf->ncmap != 0 || conf->nchist != 0) {//If there are also PNG files to delete.
				sprintf(cmd, "rm %s*.png %s", stripfname, tikzfname);
			}
			else {
				sprintf(cmd, "rm %s", tikzfname);
			}
			printf("[INFO] Removing temporary files with '%s' (pass -k to keep)...\n", cmd);
			system(cmd);
		}
		del_config(conf); //Finally, deletes all the simulation data.
	}
	return 0;  //or exit(EXIT_SUCCESS);
}
