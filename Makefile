# @author David GASPARD <dgaspard@ulb.ac.be>
# @date Created on 2020-07-26 at 12:00 CEST
# @copyright Copyright (C) 2022  David Gaspard
# @license This program is free software; it can be redistributed and/or modified under
# the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
# This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# The license file, titled "LICENSE", can be found in the root directory of this project.
# Usage:
# make all     #Compile all the binaries (msmain+msplot), but not the test programs
# make msmain  #Compile the main computation program
# make msplot  #Compile the plot program
# make test    #Compile all the test programs
# make clean   #Remove all the binaries (msmain+msplot) and objects, but not the test programs

#SOURCE FILE STRUCTURE:
SRCCOREDIR := ./src/core
SRCMAINDIR := ./src/main
SRCTESTDIR := ./src/test
BINCOREDIR := ./bin/core
BINMAINDIR := ./bin/main
BINTESTDIR := ./bin/test

SRCLST  := $(wildcard $(SRCCOREDIR)/*.c)
BINLST  := $(patsubst $(SRCCOREDIR)/%.c,$(BINCOREDIR)/%.o,$(SRCLST))

#COMPILER OPTIONS:
CC      := gcc
CFLAG   := -W -Wall -Wextra -fms-extensions -fopenmp
INCL    := -I$(SRCCOREDIR) -I$(SRCTESTDIR)
LDFLAG  := -lm -llapack

######################
# MAIN MAKE FUNCTIONS
######################
all: directories msmain msplot

directories: $(BINCOREDIR) $(BINMAINDIR) $(BINTESTDIR) 

$(BINCOREDIR) $(BINMAINDIR) $(BINTESTDIR):
	mkdir -p $@

msmain: $(BINLST) $(BINMAINDIR)/ms_main.o
	$(CC) $(CFLAG) $^ $(LDFLAG) -o $@

msplot: $(BINLST) $(BINMAINDIR)/ms_plot.o 
	$(CC) $(CFLAG) $^ $(LDFLAG) -o $@

clmain: $(BINCOREDIR)/common_util.o $(BINMAINDIR)/cl_main.o
	$(CC) $(CFLAG) $^ $(LDFLAG) -o $@

bin/%.o: src/%.c
	$(CC) $(CFLAG) $(INCL) -c $< -o $@

#################
# TEST FUNCTIONS 
#################
test: directories test_green_bessel_old test_green_bessel test_nr_random test_common_util test_vector_util test_root_finder test_medium test_ms_matrix test_complex_map test_config_io test_fitting_model test_classical_diffusion test_square_well_resonances test_eigenstate_util test_plot_functions

test_green_bessel_old: $(BINCOREDIR)/green_bessel.o $(BINTESTDIR)/test_green_bessel_old.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_green_bessel: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/complex_vector_util.o $(BINCOREDIR)/green_bessel.o $(BINTESTDIR)/assertion.o $(BINTESTDIR)/test_green_bessel.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_incomplete_gamma: $(BINCOREDIR)/incomplete_gamma.o $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/complex_vector_util.o $(BINTESTDIR)/assertion.o $(BINTESTDIR)/test_incomplete_gamma.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_nr_random: $(BINCOREDIR)/nru_random.o $(BINTESTDIR)/test_nr_random.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_common_util: $(BINCOREDIR)/common_util.o $(BINTESTDIR)/test_common_util.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_vector_util: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/complex_vector_util.o $(BINTESTDIR)/assertion.o $(BINTESTDIR)/test_vector_util.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_root_finder: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/complex_vector_util.o $(BINCOREDIR)/domain_util.o $(BINCOREDIR)/root_finder.o $(BINTESTDIR)/assertion.o $(BINTESTDIR)/test_root_finder.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_medium: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/nru_random.o $(BINCOREDIR)/green_bessel.o $(BINCOREDIR)/medium_util.o $(BINTESTDIR)/test_medium.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_ms_matrix: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/real_vector_util.o $(BINCOREDIR)/complex_vector_util.o $(BINCOREDIR)/nru_random.o $(BINCOREDIR)/medium_util.o $(BINCOREDIR)/green_bessel.o $(BINCOREDIR)/ms_matrix_util.o $(BINTESTDIR)/assertion.o $(BINTESTDIR)/test_ms_matrix.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_complex_map: $(BINLST) $(BINTESTDIR)/test_complex_map.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_config_io: $(BINLST) $(BINTESTDIR)/test_config_io.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_fitting_model: $(BINCOREDIR)/fitting_model_util.o $(BINTESTDIR)/test_fitting_model.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_classical_diffusion: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/classical_diffusion_util.o $(BINTESTDIR)/test_classical_diffusion.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_square_well_resonances: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/domain_util.o $(BINCOREDIR)/nru_random.o $(BINCOREDIR)/medium_util.o $(BINCOREDIR)/green_bessel.o $(BINCOREDIR)/square_well_resonances.o $(BINTESTDIR)/test_square_well_resonances.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_eigenstate_util: $(BINCOREDIR)/common_util.o $(BINCOREDIR)/nru_random.o $(BINCOREDIR)/green_bessel.o $(BINCOREDIR)/medium_util.o $(BINCOREDIR)/eigenstate_util.o  $(BINTESTDIR)/test_eigenstate_util.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

test_plot_functions: $(BINLST) $(BINTESTDIR)/test_plot_functions.o
	$(CC) $(CFLAG) $(INCL) $^ $(LDFLAG) -o $@

#################
# CLEAN ALL
#################
clean:
	rm -rfv bin/ msmain msplot
