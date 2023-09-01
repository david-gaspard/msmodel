/****
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created on 2020-08-03 at 13:03 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Core functions of the multiple scattering (MS) solver program.
 * Builds a rectangular medium consisting of a given number of point-like atoms placed at random.
 * Constructs the associated multiple scattering M matrix in column major representation, and returns its properties.
 ***/
#include <stdlib.h>                   /* Standard Library for Memory Allocation */
#include <stdio.h>                    /* Standard Library for Input and Output */
#include <math.h>                     /* Standard Library for Mathematical Functions */
#include <time.h>                     /* Standard Library for Time Measurement (to be replaced with OpenMP) */
#include <omp.h>                      /* Import the OpenMP Library for Multi-Processing */
#include "common_util.h"              /* Import General Purpose Functions such as the Progress Bar */
#include "real_vector_util.h"         /* Import the Real Vector Manipulation Utilities */
#include "complex_vector_util.h"      /* Import the Complex Vector Manipulation Utilities */
#include "domain_util.h"              /* Import the Complex Domain Management Utilities */
#include "green_bessel.h"             /* Import the Bessel Functions and Green Functions */
#include "ms_matrix_util.h"           /* Import the Complex Library and the Multi-Scattering Matrix Utilities */
#include "root_finder.h"              /* Import the Root-Finding Utilities to search for individual eigenstates */
#include "medium_util.h"              /* Import the Random Medium Management Utilities */
#include "config_io.h"                /* Import the Input Parsing Functions */
#include "complex_map_util.h"         /* Import the Complex Map Management Utilities */
#include "complex_cut_util.h"         /* Import the Complex Cut Management Utilities */
#include "log_cut_util.h"             /* Import the Log Cut Management Utilities */
#include "chistogram_util.h"          /* Import the Complex Histogram Utilities */
#include "eigenstate_util.h"          /* Import the Eigenstates Utilities */
#include "eigenstate_list_util.h"     /* Import the Eigenstates List Utilities */
#include "random_function_util.h"     /* Import the Random Function Utilities */
#include "diff_cross_section_util.h"  /* Import the Differential Cross Section Utilities */
#include "imag_mu_histogram_util.h"   /* Import the Utilities of the Histogram of the Imaginary Part of Mu */
#include "incident_wave_util.h"       /* Import the Incident Wave Utilities */
#include "wave_function_util.h"       /* Import the Wave Function Utilities */
#include "wave_function_cut_util.h"   /* Import the Wave Function Cut Utilities */
/**
 * Defines some macros:
 */
#define PI       3.1415926535897932385
/**
 * LAPACK's routine to solve a complex linear system using the LU decomposition method:
 */
void zgesv_(int* n, int* nrhs, dcomplex* a, int* lda, int* ipiv, dcomplex* b, int* ldb, int* info);

/*********************************************
 * COMPUTATION OF COMPLEX MAP OF THE K-PLANE
 ********************************************/
/**
 * Computes the trace-log function Tr[ln(M)]/N averaged over "ngeom" geometries in a portion of the complex k-plane, and writes the results to the complex map "cmap".
 */
void run_tracelog_serial(Medium* med, ComplexData* cdata) {
	char msg[50];    //Message to display while running.
	//char outfname[50];  //Output filename.
	dcomplex* matrix = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex)); //List of the MS matrix elements in column major order (omp-private).
	dcomplex k, trlog; //Current complex wavenumber of trace-log value (both omp-private).
	uint64_t ngeom = cdata->nseed; //The number of geometries is given by the "nseed" field found in the input file and stored in med->seed.
	uint64_t s; //Current seed, index of the current atomic geometry (omp-private).
	int j = 0;  //Current number of done jobs (omp-shared).
	int jtot = ngeom*cdata->ntot;  //Compute the total number of jobs.
	int l; //Index of the point in the complex k-plane (omp-private).
	clock_t begin = clock(); //Captures the beginning time.
	for (s = 1; s <= ngeom; s++) {//Loop over the atomic geometries.
		sprintf(msg, "TraceLog, geom: %3lu/%3lu", s, ngeom);  //Updates the message.
		fill_medium(med, s);  //Fills the medium with the random positions specified by the seed "s".
		//sprintf(fname, "%s_atoms_s%lu.tab", uuid, seed);
		//print_medium(med, fname);  //Prints the positions of the atoms to a file.
		for (l = 0; l < cdata->ntot; l++) {//Loop over the points of the complex k-plane (horizontally line by line). Typically heavy, should be parallelized with OpenMP.
			k = cdata->get_arg(cdata, l);
			build_ms_matrix(med, k, matrix); //Builds the matrix with the current value of k.
			trlog = tracelog_lu(med->n, matrix);
			cdata->data[l] = ((s-1)*cdata->data[l] + trlog)/s; //Updates the average value.
			j++; //Number of done jobs since the beginning.
			if (j%10 == 0) {//Skip a certain number of samples before refreshing the progress bar.
				print_progress(j, jtot, (double)(clock() - begin)/CLOCKS_PER_SEC, msg); //WARN: Use omp_get_wtime() with OpenMP instead.
			}
		}
		//sprintf(outfname, "%s_s%lu.dat", uuid, s);
		//cmap->realtime = (double)(clock() - begin)/CLOCKS_PER_SEC;  //Saves the total computation time. WARN: Use omp_get_wtime() with OpenMP instead.
		//save_data(cmap, med, outfname); //Saves the current map to file after each completed geometry.
	}
	free(matrix);
	cdata->realtime = (double)(clock() - begin)/CLOCKS_PER_SEC;  //Saves the total computation time. WARN: Use omp_get_wtime() with OpenMP instead.
	print_end(jtot, cdata->realtime, msg);
	//sprintf(outfname, "%s.dat", uuid);
	//save_data(cmap, med, outfname); //Saves the current map to file after each completed geometry.
}
/**
 * Computes the trace-log function Tr[ln(M)]/N averaged over "ngeom" geometries in a portion
 * of the complex k-plane, and writes the results to the complex data structure "cdata".
 * This version of map_tracelog uses OpenMP to parallelize the computation of each complex map.
 */
void run_tracelog_omp(Medium* med, ComplexData* cdata, int nthread) {
	
	uint64_t ngeom = cdata->nseed;     //The number of geometries is given by the "nseed" field from the input file and stored in med->seed (om-shared constant).
	int j = 0;                         //Current number of done jobs (omp-shared variable).
	int jtot = cdata->ntot*cdata->nseed; //Compute the total number of jobs, i.e., number of pixels times number of geometries (omp-shared constant).
	int dj = jtot/2000;                //Step size to get 0.05% precision of the progress bar (omp-shared constant).
	if (dj == 0) dj = 1;               //The step must not be zero.
	char msg[60];                      //Message to display while running (omp-shared variable).
	double begin = omp_get_wtime();    //Captures the beginning time (omp-shared constant).
	dcomplex* matrix;
	
	#pragma omp parallel shared(cdata, med, j, jtot, dj, begin, msg) private(matrix) num_threads(nthread)
	{   
		int l;           //Private index of position in the k-plane.
		uint64_t s;      //Private index of geometry.
		dcomplex k;      //Private position in the k-plane.
		dcomplex trlog;  //Private result of the trace-log function.
		matrix = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex));  //Private MS matrix.
		
		for (s = 1; s <= ngeom; s++) {//Loop on atomic geometries. Same for each thread.
			#pragma omp single
			{
				fill_medium(med, s);  //Fill the medium only once, by one thread, as it is omp-shared. Also more safe since using static random generator.
				sprintf(msg, "TraceLog, %d thr, %4lu/%4lu seed", nthread, s, ngeom);  //Updates the message.
			}
			//Cyclic/static schedule allows better load balancing in case of load increase in certain region of the k-plane. //schedule(static,1)
			#pragma omp for schedule(static,1)
			for (l = 0; l < cdata->ntot; l++) {//Loop on regions of the k-plane (parallelized).
				k = cdata->get_arg(cdata, l);
				build_ms_matrix(med, k, matrix);  //Builds the matrix with the current value of k.
				trlog = tracelog_lu(med->n, matrix);
				cdata->data[l] = ((s-1)*cdata->data[l] + trlog)/s;  //Updates the average value. No collision between threads.
				if (omp_get_thread_num() == 0) {//Print progress bar with only the main thread.
					j += nthread;  //As j is omp-shared, incrementing it out of a critical section is not thread-safe.
					if (j%dj == 0) {//This condition is met at most once per thread.
						print_progress(j, jtot, omp_get_wtime()-begin, msg); 
					}
				}
			}
			//Implicit barrier at the end of omp-for, it ensures synchronization on "med", but the process is slightly slowed down (try num_threads = num_cores).
		}
		free(matrix); //De-allocate private matrix.
	}
	cdata->realtime = omp_get_wtime()-begin;  //Saves the total computation time.
	print_end(jtot, cdata->realtime, msg);
}

/*******************************************************
 * COMPUTATION OF EIGENVALUES HISTOGRAM IN THE MU-PLANE
 ******************************************************/
/**
 * Computes the histogram of the eigenvalues of the multi-scattering matrix in the complex mu-plane for a given value of complex wave number chist->k.
 * @deprecated Old version which is not up-to-date anymore.
 */
void run_chistogram_serial(Medium* med, Chistogram* chist) {
	uint64_t s, ngeom = chist->nseed;   //Number of geometries.
	int nmu = chist->nseed*med->n;      //Total number of computed eigenvalues.
	dcomplex* matrix = (dcomplex*)malloc(med->n*med->n*sizeof(dcomplex));  //Private MS matrix (not packed because eigenvalue not supported in LAPACK).
	dcomplex* mu = (dcomplex*)calloc(nmu, sizeof(dcomplex)); //Shared list of all mu eigenvalues computed (should be <10^7).
	char msg[50] = "MuHistogram, serial";    //Info message.
	double begin = omp_get_wtime(); //Captures the beginning time (omp-shared constant).
	for (s = 1; s <= ngeom; s++) {//Loop on atomic geometries. Same for each thread.
		fill_medium(med, s);                  //Fill the medium with the given seed. Note that this is thread unsafe.
		build_ms_matrix(med, chist->k, matrix);  //Builds the (unpacked) square MS matrix with the current value of k.
		eigvals(med->n, matrix, mu + (s-1)*med->n);  //Computes the eigenvalues using LAPACK.
		print_progress(s, ngeom, omp_get_wtime()-begin, msg);
	}
	print_end(ngeom, omp_get_wtime()-begin, msg);
	int lost = fill_chist(chist, nmu, mu);   //Set up the histogram domain, and appends the eigenvalues to the histogram, after having set up the domain.
	printf("[INFO] Lost eigenvalues: %d/%d (%.2f%%)\n", lost, nmu, 100.0*lost/nmu);
	free(mu);      //De-allocate list of eigenvalues.
	free(matrix);  //De-allocate private matrix.
}
/**
 * Computes the histogram of the eigenvalues of the multi-scattering matrix in the complex mu-plane for a given value of complex wave number chist->k.
 * Use OpenMP to parallelize the different seeds.
 */
void run_chistogram_omp(Medium* med, Chistogram* chist, int nthread) {
	
	uint64_t ngeom = chist->nseed;   //Number of geometries.
	int nmu = chist->nseed*med->n;   //Total number of computed eigenvalues.
	dcomplex* mu = (dcomplex*)calloc(nmu, sizeof(dcomplex)); //OMP-shared list of all mu eigenvalues computed (should be <10^7).
	int j = 0;                          //Current number of done jobs (OMP-shared variable), to be protected against race conditions.
	char msg[50];
	sprintf(msg, "MuHistogram, %d thr", nthread);  //Info message.
	double begin = omp_get_wtime();     //Captures the beginning time (OMP-shared constant).
	Medium* locmed;
	dcomplex* matrix;
	
	#pragma omp parallel shared(chist, med, mu, nmu, j, begin, msg) private(locmed, matrix) num_threads(nthread)
	{   
		matrix = (dcomplex*)calloc(med->n*med->n, sizeof(dcomplex));  //OMP-Private (unpacked) MS matrix.
		locmed = (Medium*)calloc(1, sizeof(Medium));  //OMP-Private medium with different seed on each thread.
		init_copy_medium(med, locmed);  //Initialize a deep copy of the medium to get separate media in each thread.
		uint64_t s;  //OMP-Private index of geometry or seed.
		
		#pragma omp for schedule(static,1)
		for (s = 1; s <= ngeom; s++) {//Parallelized loop on medium geometries. Same for each thread.
			fill_medium(locmed, s);  //Fill the medium with the given seed.
			if (chist->mtype == msmatrix) {
				build_ms_matrix(locmed, chist->k, matrix);  //Builds the matrix with the current value of k.
			}
			else {//Normalized MS matrix
				build_normalized_matrix(locmed, chist->k, matrix); //Builds the normalized MS matrix which is independent from the single-atom scattering model.
			}
			eigvals(locmed->n, matrix, mu + (s-1)*locmed->n);  //Computes the eigenvalues using LAPACK on separate slots of "mu".
			#pragma omp critical(progress3)
			{
				j++;  //As j is a OMP-shared variable, incrementing it out of a critical section is not thread-safe.
				if (j%nthread == 0) {//This condition is met at most once per thread.
					print_progress(j, ngeom, omp_get_wtime()-begin, msg);
				}
			}
		}
		del_medium(locmed); //Delete the content of the OMP-private medium.
		free(locmed); //De-allocate the OMP-private medium.
		free(matrix); //De-allocate private matrix.
	}
	int lost = fill_chist(chist, nmu, mu);   //Set up the histogram domain, and appends the eigenvalues to the histogram, after having set up the domain.
	chist->realtime = omp_get_wtime()-begin; //Saves the total computation time.
	print_end(ngeom, chist->realtime, msg);
	printf("[INFO] Lost eigenvalues: %d/%d (%.2f%%)\n", lost, nmu, 100.0*lost/nmu);
	free(mu);  //De-allocate list of eigenvalues.
}

/**
 * Computes the histogram of the imaginary parts of the eigenvalues of the multi-scattering matrix for a given value of wave number "k".
 * Use OpenMP to parallelize the different seeds.
 */
void run_imag_mu_histogram_omp(Medium* med, ImagMuHistogram* ihist, int nthread) {
	
	double begin = omp_get_wtime();  //Captures the beginning time (OMP-shared constant).
	uint64_t ngeom = ihist->nseed;   //Number of geometries (OMP-Shared constant).
	int j = 0;     //Current number of done jobs (OMP-shared variable), to be protected against race conditions.
	setup_domain_ihist(ihist, med); //Use the seed s=1 to initialize the histogram bounds.
	int lost = 0;  //Number of lost eigenvalues (OMP-Shared variable), to be protected against race conditions.
	ihist->avgx = creal(free_green_imag(med->d, ihist->k, 0.));  //Compute the exact average of Im(mu). This is the exact average for any possible random realization of the medium.
	char msg[60];
	sprintf(msg, "Im(Mu) Histogram, %d thr", nthread);  //Info message.
	Medium* locmed;
	dcomplex *matrix, *mu;
	
	#pragma omp parallel shared(ihist, j, lost, begin, msg) private(matrix, mu, locmed) num_threads(nthread)
	{   
		matrix = (dcomplex*)calloc(med->n*med->n, sizeof(dcomplex));  //OMP-Private (unpacked) MS matrix.
		mu = (dcomplex*)calloc(med->n, sizeof(dcomplex)); //OMP-Private list of eigenvalues.
		locmed = (Medium*)calloc(1, sizeof(Medium));  //OMP-Private medium with different seed on each thread.
		init_copy_medium(med, locmed);  //Initialize a OMP-private deep copy of the medium to get separate media in each thread.
		double varx_true = 0., varx_expc = 0.; //True and expected variances of Im(mu) (OMP-Private).
		uint64_t s;  //OMP-Private index of geometry or seed.
		int i; //OMP-Private loop index.
		
		#pragma omp for schedule(static,1)
		for (s = 1; s <= ngeom; s++) {//Parallelized loop on medium geometries. Same for each thread.
			fill_medium(locmed, s);   //Fill the medium with the given seed.
			build_ms_matrix(locmed, ihist->k, matrix);  //Builds the matrix with the current value of k.
			eigvals(locmed->n, matrix, mu);  //Computes the eigenvalues using LAPACK on private slots "mu". Most time-consuming operation.
			varx_expc += 0.25*tracegtwo(locmed->n, matrix); //Compute Tr[G⁺G]/N from M to estimate the variance of Im(mu).
			for (i = 0; i < locmed->n; i++) {//Loop on the eigenvalues to compute the true variance of Im(mu). This loop does not need parallelization.
				varx_true += (cimag(mu[i]) - ihist->avgx)*(cimag(mu[i]) - ihist->avgx);
			}
			#pragma omp critical(merge_eigvals)
			{
				lost += append_list_ihist(ihist, locmed->n, mu); //Each thread should append its eigenvalues separately (one thread at a time).
				j++;  //Update the current number of done jobs.
				if (j%nthread == 0) //This condition is met at most once per thread.
					print_progress(j, ngeom, omp_get_wtime()-begin, msg);
			}
		}
		#pragma omp critical(merge_moments)
		{
			ihist->varx_true += varx_true; //Merges/reduces the variances.
			ihist->varx_expc += varx_expc;
		}
		del_medium(locmed); //Delete the content of the private medium.
		free(locmed); //De-allocate the private medium.
		free(matrix); //De-allocate private matrix.
		free(mu); //De-allocate the private eigenvalues.
	}
	normalize_ihist(ihist); //Finally, normalizes the histogram to get a probability density.
	int nmu = ihist->nseed*med->n; //Total number of eigenvalues.
	ihist->varx_true /= nmu;  //Normalize the variances.
	ihist->varx_expc /= ihist->nseed;
	ihist->realtime = omp_get_wtime()-begin; //Saves the total computation time.
	print_end(ngeom, ihist->realtime, msg);
	printf("[INFO] Lost eigenvalues: %d/%d (%.2f%%)\n", lost, nmu, 100.0*lost/nmu);
}

/*********************************************
 * COMPUTATION OF EIGENSTATES IN THE K PLANE
 ********************************************/
/**
 * Initializes and computes the "nroot" eigenstates at the given positions "kroot" using the inverse power iteration with the prescribed parameters (maxit/toler/verb).
 * @param "eigls" Eigenstate list structure where the eigenvectors data are saved.
 * @param "nroot" Number of found roots in the complex k plane.
 * @param "kroot" Positions of the found roots in the complex k plane.
 * @param "maxit" Maximum number of iterations of the inverse power iteration (typically 50).
 * @param "toler" Tolerance of the relative distance between two successive vectors obtained by inverse iterations (typically 1e-10).
 * @param "verb" Verbosity level of the inverse iteration function (0=quiet, 1=verbose).
 */
void compute_eigenstates(EigenstateList* eigls, int nroot, dcomplex* kroot, int maxit, double toler, int verb) {
	int i, n = eigls->med->n;
	dcomplex mu, *matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	alloc_eigenstate_list(eigls, nroot); //Allocate space for all the complex resonance roots that have been found.
	for (i = 0; i < eigls->nstate; i++) {//Loop on the roots to compute the associated eigenvector.
		printf("\r[EXEC] Computing state %3d/%3d, kroot=%20.14g%+20.14gi ... ", (i+1), eigls->nstate, creal(kroot[i]), cimag(kroot[i]));
		if (verb >= 1) printf("\n");
		fflush(stdout);
		eigls->eigtab[i]->kroot = kroot[i];  //Save the position of the complex resonance root "kroot".
		build_ms_matrix(eigls->med, kroot[i], matrix); //Builds the MS matrix at "kroot" in prevision of the inverse iteration.
		constant_unit_cvector(n, eigls->eigtab[i]->eigv); //Initialized a normalized complex vector.
		inverse_iteration(n, matrix, eigls->eigtab[i]->eigv, &mu, maxit, toler, verb); //Execute a few inverse iterations to compute the eigenvector/eigenvalue.
		eigls->eigtab[i]->mu = mu; //Saves the corresponding eigenvalue.
	}
	printf("Done\n");
	check_eigenstate_list(eigls);  //Check that all eigenstates are valid.
	free(matrix);
} 

/**
 * Computes the eigenstates planned in the given eigenstate list structure.
 * First find distinct complex roots of det(M(k))=0 using the Maehly-Aberth-Ehrlich iteration,
 * then determine the corresponding eigenstates using a few steps of inverse power iteration.
 */
void find_eigenstates_serial_1(EigenstateList* eigls, int* nroot, dcomplex* kroot) {
	dcomplex dlogf(dcomplex k) {//Defines the logarithmic derivative, f'(k)/f(k), of the function to solve, i.e., f(k)=det(M(k)).
		return dk_log_det_ms(eigls->med, k);
	};
	double begin = omp_get_wtime();  //Captures the beginning time (OMP-shared constant).
	int nexpsup = 4; //Default number of parameters involved in the exponential suppression.
	solve_aberth(dlogf, eigls->dom, nroot, kroot, nexpsup, eigls->maxit, eigls->toler, eigls->verb); //Finds several distinct roots of det(M(k))=0 using Aberth method.
	eigls->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
}

/**
 * Computes the eigenstates planned in the given eigenstate list structure.
 * First find distinct complex roots of det(M(k))=0 using the secant iteration on the minimum eigenvalue of M(k),
 * then determine the corresponding eigenstates using a few steps of inverse power iteration.
 */
void find_eigenstates_serial_2(EigenstateList* eigls, int* nroot, dcomplex* kroot) {
	int* conv = (int*)calloc(*nroot, sizeof(dcomplex)); //Convergence flag for each root (0=failed, 1=converged).
	int n = eigls->med->n; //Number of atoms and size of the MS matrix.
	dcomplex* matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	dcomplex* vec = (dcomplex*)calloc(n, sizeof(dcomplex)); //Local eigenvector.
	constant_unit_cvector(n, vec); //Normalize the complex vector used by f(k) only once to reduce as many spurious zero eigenvalues as possible.
	int ivit_maxit = 50;  //Maximum number of iterations of the inverse power iteration.
	double ivit_toler = 1e-10;  //Tolerance of the relative distance between two successive vectors obtained by inverse iterations.
	int ivit_verb = 0;  //Verbosity level of the inverse iteration function (0=quiet, 1=verbose).
	dcomplex mineigval(dcomplex k) {//Defines the minimum eigenvalue function f(k), which is actually a piecewise complex function.
		dcomplex mu; //Smallest eigenvalue.
		build_ms_matrix(eigls->med, k, matrix); //Builds the MS matrix at "k" in anticipation of the inverse iteration.
		inverse_iteration(n, matrix, vec, &mu, ivit_maxit, ivit_toler, ivit_verb); //Execute a few inverse iterations to compute the eigenvector/eigenvalue.
		return mu;
	};
	double begin = omp_get_wtime(); //Captures the beginning time (OMP-shared constant).
	char msg[30];
	int i;
	for (i = 0; i < *nroot; i++) {//Loop on the target roots.
		sprintf(msg, "Root %4d/%4d", (i+1), *nroot);
		conv[i] = find_root_muller(mineigval, eigls->dom, &kroot[i], eigls->maxit, eigls->toler, eigls->verb); //Find te roots of f(k) using the secant method.
		conv[i] *= (cabs(mineigval(kroot[i])) < sqrt(eigls->toler)); //In addition, the minimum eigenvalue must be small enough to declare convergence.
		print_progress((i+1), *nroot, omp_get_wtime()-begin, msg);
		if (eigls->verb >= 1) printf("\n");
	}
	print_end(*nroot, omp_get_wtime()-begin, msg);
	int nfail = filter_cvector(nroot, kroot, conv); //Removes the non-converged roots.
	double dupl_toler = 1e-14;  //Relative tolerance between different roots (when eliminating duplicates).
	int ndupl = uniq_cvector(nroot, kroot, dupl_toler); //Removes duplicate roots with relatively large scope radius.
	printf("[INFO] Found %d/%d valid roots, %d duplicates, %d failed.\n", *nroot, eigls->ntarget, ndupl, nfail);
	eigls->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	free(conv);
	free(matrix);
	free(vec);
}

/**
 * Computes the eigenstates planned in the given eigenstate list structure.
 * First find distinct complex roots of det(M(k))=0 using the secant iteration on the function 1/Tr[M(k)^-1],
 * then determine the corresponding eigenstates using a few steps of inverse power iteration.
 */
void find_eigenstates_serial_3(EigenstateList* eigls, int* nroot, dcomplex* kroot) {
	int* conv = (int*)calloc(*nroot, sizeof(dcomplex)); //Convergence flag for each root (0=failed, 1=converged).
	int n = eigls->med->n; //Number of atoms and size of the MS matrix.
	dcomplex* matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));
	dcomplex invtrinv(dcomplex k) {//Defines the function 1/Tr(M^-1) for which the roots are sought.
		build_ms_matrix(eigls->med, k, matrix); //Builds the MS matrix at "k".
		return 1./traceinv(n, matrix);
	}
	double begin = omp_get_wtime(); //Captures the beginning time (OMP-shared constant).
	char msg[30];
	int i;
	for (i = 0; i < *nroot; i++) {//Loop on the target roots.
		sprintf(msg, "Root %4d/%4d", (i+1), *nroot);
		conv[i] = find_root_muller(invtrinv, eigls->dom, &kroot[i], eigls->maxit, eigls->toler, eigls->verb); //Find te roots of f(k) using the secant method.
		print_progress((i+1), *nroot, omp_get_wtime()-begin, msg);
		if (eigls->verb >= 1) printf("\n");
	}
	print_end(*nroot, omp_get_wtime()-begin, msg);
	int nfail = filter_cvector(nroot, kroot, conv); //Removes the non-converged roots.
	double dupl_toler = 1e-14;  //Relative tolerance between different roots (when eliminating duplicates).
	int ndupl = uniq_cvector(nroot, kroot, dupl_toler); //Removes duplicate roots with relatively large scope radius.
	printf("[INFO] Found %d/%d valid roots, %d duplicates, %d failed.\n", *nroot, eigls->ntarget, ndupl, nfail);
	eigls->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	free(conv);
	free(matrix);
}

/**
 * Computes the eigenstates planned in the given eigenstate list structure.
 * Select among the available method to find the roots of det(M(k))=0 depending on the user choice.
 */
void run_eigenstate_list_serial(EigenstateList* eigls) {
	printf("[INFO] Computing %d roots in region x=%g:%g y=%g:%g, with %d iterations of method '%s'.\n",
		eigls->ntarget, eigls->dom->xmin, eigls->dom->xmax, eigls->dom->ymin, eigls->dom->ymax, eigls->maxit, eigmethod_to_string(eigls->method));
	int nroot = eigls->ntarget; //Initial number of roots.
	dcomplex* kroot = (dcomplex*)calloc(nroot, sizeof(dcomplex)); //Allocate space for the complex roots of det(M(k))=0.
	set_halton_points(eigls->dom, nroot, kroot); //Initialize the roots with quasi-random points in the complex domain.
	double opposite_imag(dcomplex z) { return -cimag(z); }; //Opposite imaginary part to sort the roots from closest to real axis to far away.
	sort_cvector(nroot, kroot, opposite_imag); //Pre-sort the roots to better exploit the eigenstate computed form the previous iteration (assuming they are close to each other).
	if (eigls->verb >= 2) {
		print_cvector(nroot, kroot, "Root guesses"); fflush(stdout);
	}
	switch (eigls->method) {
		case determinant:
			find_eigenstates_serial_1(eigls, &nroot, kroot);
			break;
		case mineigval:
			find_eigenstates_serial_2(eigls, &nroot, kroot);
			break;
		case invtraceinv:
			find_eigenstates_serial_3(eigls, &nroot, kroot);
			break;
	}
	sort_cvector(nroot, kroot, opposite_imag); //Sort the roots according to their imaginary part.
	int ivit_maxit = 50;       //Maximum number of iterations of the inverse power iteration.
	double ivit_toler = 1e-10; //Tolerance of the relative distance between two successive vectors obtained by inverse iterations.
	int ivit_verb = 0;         //Verbosity level of the inverse iteration function (0=quiet, 1=verbose).
	compute_eigenstates(eigls, nroot, kroot, ivit_maxit, ivit_toler, ivit_verb);
	free(kroot);
}

/*********************************************
 * COMPUTATION OF THE RANDOM REAL FUNCTIONS
 ********************************************/
/**
 * Computes in serial all the samples of the random function prescribed by "rfun->type", either the total cross section or the density of states.
 * The initial particle wave function is assumed to be a plane wave oriented in the main (1x) direction.
 * The medium seed is incremented after each sample so that an atomic geometry is never encountered twice. In this way, abnormal samples do not persist.
 */
void run_random_function_serial(Medium* med, RandomFunction* rfun) {
	double begin = omp_get_wtime(); //Initial time (OMP-Shared constant).
	int jtot = rfun->ntot;  //Total number of computational jobs (OMP-Shared constant).
	char msg[30]; //Info message (OMP-Shared constant).
	sprintf(msg, "Rnd function (%s)", random_function_type_to_string(rfun->type));
	int j, i; //j=Current job index (OMP-Private). i=Current bin index (OMP-Private).
	double k; //Current real value of the wave number (OMP-Private).
	
	for (j = 0; j < jtot; j++) {//Loop on the computation jobs.
		i = j/rfun->nseed; //Current bin index.
		k = rfun->xmin + i*(rfun->xmax - rfun->xmin)/(rfun->nbin-1);
		fill_medium(med, j); //Fill the medium with the given seed.
		if (rfun->type == crsec) {//Total cross section is computed.
			rfun->data[j] = total_cross_section(med, k);
		}
		else {//Density of states per unit of "k^2" is computed.
			rfun->data[j] = -cimag(dk_log_det_ms(med, k))/(2*k*PI);
		}
		if (j%10 == 0) {
			print_progress(j+1, jtot, omp_get_wtime()-begin, msg); 
		}
	}
	print_end(jtot, omp_get_wtime()-begin, msg);
	sort_random_function(rfun); //Sorts the samples in every bin, just once.
	rfun->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
}

/**
 * Computes in parallel (with OpenMP) all the samples of the random function prescribed by "rfun->type", either the total cross section or the density of states.
 * The initial particle wave function is assumed to be a plane wave oriented in the main (1x) direction.
 */
void run_random_function_omp(Medium* med, RandomFunction* rfun, int nthread) {
	double begin = omp_get_wtime();  //Initial time (OMP-Shared constant).
	uint64_t ngeom = rfun->nseed;    //Total number of different atomic geometries (seeds).
	int j, jtot = rfun->ntot;        //Total number of computational jobs (OMP-Shared constant).
	int dj = jtot/1000;              //Step size to get 0.1% precision of the progress bar (OMP-Shared constant).
	if (dj < nthread) dj = nthread;  //The step must not be smaller than the number of threads.
	char msg[80]; //Info message (OMP-Shared constant).
	
	#pragma omp parallel shared(rfun, med, begin, j, jtot, dj, msg) num_threads(nthread)
	{
		int i;         //Private index of position in the k-plane.
		uint64_t s, g; //Private indices of geometry.
		double k;      //Current real value of the wave number (OMP-Private).
		
		for (g = 0; g < ngeom; g++) {//Loop on atomic geometries. Same for each thread.
			#pragma omp single
			{
				s = g+1; //Current seed.
				fill_medium(med, s);  //Fill the medium only once, by one thread, as it is omp-shared. Also more safe since using static random generator.
				sprintf(msg, "Rnd fct (%s), %d thr, %4lu/%4lu seed", random_function_type_to_string(rfun->type), nthread, s, ngeom);  //Updates the message.
			}
			#pragma omp for schedule(static,1)
			for (i = 0; i < rfun->nbin; i++) {//Loop on the bins.
				k = get_arg_rfun(rfun, i); //Get the wave number "k" corresponding to the given index "i".
				if (rfun->type == crsec) {//Total cross section is computed.
					rfun->data[g + i*ngeom] = total_cross_section(med, k);
				}
				else {//Density of states per unit of "k^2" is computed.
					rfun->data[g + i*ngeom] = -cimag(dk_log_det_ms(med, k))/(2*k*PI);
				}
				if (omp_get_thread_num() == 0) {//Print progress bar with only the main thread.
					j += nthread;  //As j is omp-shared, incrementing it out of a critical section is not thread-safe.
					if (j%dj == 0) {//This condition is met at most once per thread.
						print_progress(j, jtot, omp_get_wtime()-begin, msg); 
					}
				}
			}
			//Implicit barrier at the end of OMP-for.
		}
	}
	sort_random_function(rfun); //Sorts the samples in every bin, just once.
	rfun->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	print_end(jtot, omp_get_wtime()-begin, msg);
}

/**
 * Computes in parallel (with OpenMP) the differential cross section described in the "dcs" structure for many atomic geometries.
 * The initial particle wave function is assumed to be a plane wave oriented in the main (1x) direction.
 * @warning This function should be checked before use !!
 */
void run_differential_cross_section_omp(Medium* med, DiffCrossSection* dcs, int nthread) {
	if (med->d == 1)
		printf("[WARN] The differential cross section in %dD cannot be computed in the same way as in other dimensions. This may lead to strange results...\n", med->d);
	double begin = omp_get_wtime(); //Initial time (OMP-Shared constant).
	int ngeom = dcs->nseed;  //Number of atomic geometries.
	int ds = ngeom/1000;     //Step size to get 0.1% precision of the progress bar (OMP-Shared constant).
	if (ds < nthread) ds = nthread;    //The step must not be smaller than the number of threads.
	char msg[] = "Diff cross section"; //Info message (OMP-Shared constant).
	int n = med->n, nrhs = 1, *ipiv;  //nrhs=Number of right-hand sides used while solving the linear system of equations for the amplitudes, M*a = phi.
	double prefac = pow(dcs->k/(2*PI), med->d-3)/(16*PI*PI); //Prefactor used to computes the differential cross section from the scattering amplitudes.
	double* samples = (double*)calloc(ngeom*dcs->nbin, sizeof(double));  //All the samples of differential cross section (OMP-Shared variable).
	Medium* locmed; //Current medium which is different for each thread (OMP-Private).
	dcomplex *matrix, *a, *b; //matrix=Current MS matrix. a=Current vector of scattering amplitudes per atom. b=Current outgoing plane wave vector (all OMP-Private).
	
	#pragma omp parallel shared(dcs, samples, begin, msg) private(locmed, matrix, ipiv, a, b) num_threads(nthread)
	{	
		locmed = (Medium*)calloc(1, sizeof(Medium));  //Local medium with different seed on each thread (OMP-Private).
		init_copy_medium(med, locmed);  //Initialize a deep copy of the medium to get separate media in each thread.
		matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));  //Allocate the MS matrix (OMP-Private).
		ipiv = (int*)calloc(n, sizeof(int));  //Pivot index array used by LAPACK (OMP-Private).
		a = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of scattering amplitudes (OMP-Private).
		b = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of output plane wave vector (OMP-Private).
		int s, i, info; //s=Current seed, i=Current bin index, info=Info flag used by LAPACK (all OMP-Private).
		dcomplex atot; //Total scattering amplitude in a given direction (OMP-Private).
		double theta; //Current scattering angle (OMP-Private).
		
		#pragma omp for schedule(static,1)
		for (s = 0; s < ngeom; s++) {//Loop on the atomic geometries.
			fill_medium(locmed, s + dcs->iseed); //Reset the atomic geometry with seed "s".
			build_ms_matrix(locmed, dcs->k, matrix); //Builds the matrix with the given value of "k".
			build_plane_wave(locmed, dcs->k, 0., a); //Builds the plane wave oriented in the 1x direction (theta=0).
			zgesv_(&n, &nrhs, matrix, &n, ipiv, a, &n, &info); //Solve the linear system M*a = phi, so that "a" now contains the scattering amplitudes per atom.
			for (i = 0; i < dcs->nbin; i++) {//Loop on the scattering angles, one per bin.
				theta = get_arg_dcs(dcs, i); //Get the current abscissa angle in degrees.
				build_plane_wave(locmed, dcs->k, theta, b); //Builds the plane wave oriented in the output direction "th" (in degrees).
				atot = cdot_product(n, b, a); //Deduce the total scattering amplitude.
				samples[s + i*ngeom] = prefac*(creal(atot)*creal(atot) + cimag(atot)*cimag(atot)); //Computes the differential cross section.
			}
			if (s%ds == 0) //Display the progress bar at most once per thread.
				print_progress(s+1, ngeom, omp_get_wtime()-begin, msg);
		}
		del_medium(locmed);
		free(locmed); free(matrix);
		free(a); free(b); free(ipiv);
	}
	int s, i, q;
	for (s = 0; s < dcs->nc; s++) {//Loop on the first sample curves to plot them as an overview.
		for (i = 0; i < dcs->nbin; i++) {//Copy the bins of the s^th curve.
			dcs->curve[s + i*dcs->nc] = samples[s + i*ngeom];
		}
	}
	for (i = 0; i < dcs->nbin; i++) {//Loop on the bins to compute the quantiles and the mean.
		sort_vector(ngeom, samples + i*ngeom); //The samples must be sorted before computing quantiles (this may take some time).
		for (q = 0; q < dcs->nq; q++) {//Computes the quantiles for each bin.
			dcs->data[q + i*dcs->nq] = quantile(ngeom, samples + i*ngeom, ((double)q)/(dcs->nq-1));
		}
		dcs->mean[i] = mean_vector(ngeom, samples + i*ngeom); //Compute the mean curve.
	}
	free(samples);
	dcs->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	print_end(ngeom, omp_get_wtime()-begin, msg);
}

/*********************************************
 * COMPUTATION OF THE AVERAGE WAVE FUNCTION
 ********************************************/

/**
 * Computes the average wave function described in the "wfun" structure for many atomic geometries.
 */
void run_wave_function_omp(Medium* med, WaveFunction* wfun, int nthread) {
	if (med->d == 1)
		printf("[WARN] The wavefunction in %dD cannot be computed in the same way as in other dimensions. This may lead to strange results...\n", med->d);
	double begin = omp_get_wtime(); //Initial time (OMP-Shared constant).
	int s, ngeom = wfun->nseed, ntot = wfun->nx*wfun->ny; //Total numbers of configurations or pixels (MP-Shared constants).
	int i, d = med->d, n = med->n, nrhs = 1, info;
	int j = 0;  //Current number of done jobs (OMP-shared variable, to be protected against race conditions).
	int jtot = ntot*ngeom;  //Total number of jobs (OMP-shared constant).
	int dj = jtot/1000; //Step size to get 0.1% precision of the progress bar (OMP-Shared constant).
	if (dj < 1) dj = 1;  //The step must not be smaller than the number of threads.
	char msg[80];  //Info message (OMP-Shared in single).
	dcomplex* matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex)); //Current MS matrix (OMP-Shared in single).
	dcomplex* a = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of scattering amplitudes (OMP-Shared in single).
	int* ipiv = (int*)calloc(n, sizeof(int));  //Pivot index array used by LAPACK (OMP-Shared in single).
	double* pos; //Current position where the wave function is evaulated (OMP-Private).
	
	for (s = 0; s < ngeom; s++) {//Loop on the atomic geometries.
		//Fill the medium and solve the MS system only once, by a single thread, as it is OMP-shared.
		fill_medium(med, s + wfun->iseed);  //Reset the atomic geometry with the given seed.
		sprintf(msg, "Wavefun, %d thr, %4d/%4d seed", nthread, s+1, ngeom);  //Updates the message.
		for (i = 0; i < n; i++) {//Loop on the elements of the incident vector.
			a[i] = eval_incident_wave(wfun->iwave, med->d, med->pos+i*d);
		}
		build_ms_matrix(med, wfun->iwave->k, matrix); //Builds the matrix with the given value of "k".
		zgesv_(&n, &nrhs, matrix, &n, ipiv, a, &n, &info); //Solve the linear system M*a = phi, so that "a" now contains the scattering amplitudes per atom.
		
		#pragma omp parallel shared(wfun, med, a, msg, begin) private(pos) num_threads(nthread)
		{
			pos = (double*)calloc(d, sizeof(double)); //Set up the current position with zero values. Only a plane is covered, other directions are ignored...
			int l; //Index of current configuration (OMP-Private). Position index in space (OMP-Private).
			dcomplex psi; //Current total wave function value (OMP-Private).
			
			//Parallelized loop on picture elements. Cyclic/static schedule allows better load balancing.
			#pragma omp for schedule(static,1)
			for (l = 0; l < ntot; l++) {//Loop on the pixels of the picture.
				get_pos_wfun(wfun, l, pos); //Gets the current position in the 2D plane.
				psi = eval_total_wave(wfun->iwave, med, a, pos); //Evaluate the total wave function.
				wfun->function[l] += creal(psi); //Adds the real part of the wave function to the 2D array.
				wfun->density[l] += creal(psi)*creal(psi) + cimag(psi)*cimag(psi); //Adds the density, |psi|^2, to the 2D array.
				if (omp_get_thread_num() == 0) {//Print progress bar with only the main thread.
					j += nthread;  //As j is OMP-shared, it must be protected to be incremented only once per thread.
					if (j%dj == 0) //Limits the amount of refresh per unit time.
						print_progress(j, jtot, omp_get_wtime()-begin, msg);
				}
			}
			free(pos);
		}
	}
	free(matrix);
	free(a); free(ipiv);
	divide_wfun(wfun, ngeom); //Normalize the wave function to get the average.
	wfun->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	print_end(ngeom, wfun->realtime, msg);
}

/**
 * Computes a longitudinal cut of the wave function (average and fluctuations) described in the "wfcut" structure for many atomic geometries.
 */
void run_wave_function_cut_omp(Medium* med, WaveFunctionCut* wfcut, int nthread) {
	double begin = omp_get_wtime(); //Initial time (OMP-Shared constant).
	int ngeom = wfcut->nseed;  //Number of atomic geometries.
	int ds = ngeom/1000;     //Step size to get 0.1% precision of the progress bar (OMP-Shared constant).
	if (ds < nthread) ds = nthread;    //The step must not be smaller than the number of threads.
	char msg[] = "Wavefuncut"; //Info message (OMP-Shared constant).
	int n = med->n, d = med->d, nrhs = 1, *ipiv;  //nrhs=Number of right-hand sides used while solving the linear system of equations for the amplitudes, M*a = phi.
	double* pos; //Current position in space where the wave function is evaluated (OMP-Private).
	double* repsi_samples = (double*)calloc(ngeom*wfcut->nbin, sizeof(double));    //All the samples of the real part of the wave function (OMP-Shared variable).
	double* impsi_samples = (double*)calloc(ngeom*wfcut->nbin, sizeof(double));    //All the samples of the imaginary part of the wave function (OMP-Shared variable).
	double* density_samples = (double*)calloc(ngeom*wfcut->nbin, sizeof(double));  //All the samples of the square modulus of the wave function (OMP-Shared variable).
	Medium* locmed; //Current medium which is different for each thread (OMP-Private).
	dcomplex *matrix, *a; //matrix=Current MS matrix. a=Current vector of scattering amplitudes per atom. (all OMP-Private).
	
	#pragma omp parallel shared(wfcut, repsi_samples, impsi_samples, density_samples, begin, msg) private(locmed, matrix, ipiv, a, pos) num_threads(nthread)
	{	
		locmed = (Medium*)calloc(1, sizeof(Medium));  //Local medium with different seed on each thread (OMP-Private).
		init_copy_medium(med, locmed);  //Initialize a deep copy of the medium to get separate media in each thread.
		matrix = (dcomplex*)calloc(n*n, sizeof(dcomplex));  //Allocate the MS matrix (OMP-Private).
		ipiv = (int*)calloc(n, sizeof(int));  //Pivot index array used by LAPACK (OMP-Private).
		a = (dcomplex*)calloc(n, sizeof(dcomplex));  //Allocate the vector of scattering amplitudes (OMP-Private).
		pos = (double*)calloc(d, sizeof(double));  //Allocate the current position where the wave function is evaluated (OMP-Private).
		dcomplex psi;   //Value of the total wave function at some position (OMP-Private).
		int s, i, info; //s=Current seed, i=Current bin index, info=Info flag used by LAPACK (all OMP-Private).
		
		#pragma omp for schedule(static,1)
		for (s = 0; s < ngeom; s++) {//Loop on the atomic geometries.
			fill_medium(locmed, s + wfcut->iseed); //Reset the atomic geometry with seed "s".
			build_ms_matrix(locmed, wfcut->iwave->k, matrix); //Builds the matrix with the given value of "k".
			for (i = 0; i < locmed->n; i++) {//Loop on the scatterers, i.e. the elements of the incident vector.
				a[i] = eval_incident_wave(wfcut->iwave, d, locmed->pos+i*d); //Builds the incident wave.
			}
			zgesv_(&n, &nrhs, matrix, &n, ipiv, a, &n, &info); //Solve the linear system M*a = phi, so that "a" now contains the scattering amplitudes per atom.
			for (i = 0; i < wfcut->nbin; i++) {//Loop on the bins to compute the total wave function.
				get_pos_wfcut(wfcut, i, pos); //Gets the current position "pos" in space.
				psi = eval_total_wave(wfcut->iwave, locmed, a, pos); //Evaluate the total wave function.
				repsi_samples[s + i*ngeom] = creal(psi);
				impsi_samples[s + i*ngeom] = cimag(psi);
				density_samples[s + i*ngeom] = creal(psi)*creal(psi) + cimag(psi)*cimag(psi);
			}
			if (s%ds == 0) //Display the progress bar at most once per thread.
				print_progress(s+1, ngeom, omp_get_wtime()-begin, msg);
		}
		del_medium(locmed);
		free(locmed); free(matrix);
		free(a); free(ipiv);
	}
	store_data_wfcut(wfcut, repsi_samples, impsi_samples, density_samples); //Filters the data and store the result in "wfcut".
	free(repsi_samples); free(impsi_samples); free(density_samples); //Free allocated memory.
	wfcut->realtime = omp_get_wtime()-begin;  //Saves the total computation time (in seconds).
	print_end(ngeom, omp_get_wtime()-begin, msg);
}

/*****************
 * MAIN FUNCTION
 ****************/
/**
 * Main function of multiple scattering computation core file. Read the simulation parameters from the input data file "infname", and writes the output to the file "outfname".
 * The integer "nthread" is the number of threads used by OpenMP for the functions which allow multi-processing.
 */
void ms_core_main(Config* conf, const char* outfname, int nthread) {
	
	print_param_medium(conf->med);    //Display the values of the parameters of the medium to standard output for information.
	save_medium(conf->med, outfname); //Already save the medium to the file.
	
	int i;
	for (i = 0; i < conf->ncmap; i++) {//Loop on complex maps in the k-plane.
		printf("====== Map of the k-plane (%d/%d) \n", i+1, conf->ncmap);
		print_param_complexmap(conf->cmap+i, 0);
		detect_roundoff_error(conf->med, (conf->cmap+i)->ymin, 0); //Overflow detection before the construction of the MS matrix.
		run_tracelog_omp(conf->med, (ComplexData*)(conf->cmap+i), nthread); //Compute the trace-log function of the MS matrix.
		save_complexmap(conf->cmap+i, outfname); //Save the computed map in the output file.
	}
	for (i = 0; i < conf->nccut; i++) {//Loop on complex cuts in the k-plane.
		printf("====== Cut of the k-plane (%d/%d) \n", i+1, conf->nccut);
		print_param_complexcut(conf->ccut+i);
		print_cross_section(conf->med, (conf->ccut+i)->zmin); //Prints the cross section information at the first point.
		run_tracelog_omp(conf->med, (ComplexData*)(conf->ccut+i), nthread); //Compute the trace-log function of the MS matrix.
		save_complexcut(conf->ccut+i, outfname); //Save the computed cut in the output file.
	}
	for (i = 0; i < conf->nlcut; i++) {//Loop on log cuts in the k-plane.
		printf("====== Log cut of the k-plane (%d/%d) \n", i+1, conf->nlcut);
		print_param_logcut(conf->lcut+i);
		detect_roundoff_error(conf->med, (conf->lcut+i)->imkmax, 0); //Overflow detection before the construction of the MS matrix.
		run_tracelog_omp(conf->med, (ComplexData*)(conf->lcut+i), nthread); //Compute the trace-log function of the MS matrix.
		save_logcut(conf->lcut+i, outfname); //Save the computed log cut in the output file.
	}
	for (i = 0; i < conf->nchist; i++) {//Loop on complex histograms in the mu-plane.
		printf("====== Histogram of the mu-plane (%d/%d) \n", i+1, conf->nchist);
		print_param_chist(conf->chist+i);
		detect_roundoff_error(conf->med, cimag((conf->chist+i)->k), 0); //Overflow detection before the construction of the MS matrix.
		run_chistogram_omp(conf->med, conf->chist+i, nthread); //Compute the mu eigenvalue density.
		save_chistogram(conf->chist+i, outfname); //Save the computed histogram in the output file.
	}
	for (i = 0; i < conf->neigls; i++) {//Loop on eigenstates objects.
		printf("====== Eigenstate list (%d/%d) \n", i+1, conf->neigls);
		run_eigenstate_list_serial(conf->eigls+i);
		save_eigenstate_list(conf->eigls+i, outfname); //Save the computed eigenstates in the output file.
	}
	for (i = 0; i < conf->nrfun; i++) {//Loop on random function objects.
		printf("====== Random function (%d/%d) \n", i+1, conf->nrfun);
		print_param_random_function(conf->rfun+i);
		run_random_function_omp(conf->med, conf->rfun+i, nthread);
		save_random_function(conf->rfun+i, outfname); //Save the computed random function in the output file.
	}
	for (i = 0; i < conf->ndcs; i++) {//Loop on random function objects.
		printf("====== Differential Cross Section (%d/%d) \n", i+1, conf->ndcs);
		print_param_diff_cross_section(conf->dcs+i);
		run_differential_cross_section_omp(conf->med, conf->dcs+i, nthread);
		save_diff_cross_section(conf->dcs+i, outfname); //Save the computed differential cross section in the output file.
	}
	for (i = 0; i < conf->nihist; i++) {//Loop on histogram Im(mu) objects.
		printf("====== Histogram Im(mu) (%d/%d) \n", i+1, conf->nihist);
		print_param_ihist(conf->ihist+i);
		run_imag_mu_histogram_omp(conf->med, conf->ihist+i, nthread);
		save_ihist(conf->ihist+i, outfname); //Save the computed histogram in the output file.
	}
	for (i = 0; i < conf->nwfun; i++) {//Loop on histogram Im(mu) objects.
		printf("====== Wave Function (%d/%d) \n", i+1, conf->nwfun);
		print_param_wfun(conf->wfun+i);
		run_wave_function_omp(conf->med, conf->wfun+i, nthread);
		save_wfun(conf->wfun+i, outfname); //Save the computed histogram in the output file.
	}
	for (i = 0; i < conf->nwfcut; i++) {//Loop on histogram Im(mu) objects.
		printf("====== Wave Function Cut (%d/%d) \n", i+1, conf->nwfcut);
		print_param_wfcut(conf->wfcut+i);
		run_wave_function_cut_omp(conf->med, conf->wfcut+i, nthread);
		save_wfcut(conf->wfcut+i, outfname); //Save the computed histogram in the output file.
	}
	printf("====== Data saved to '%s' \n", outfname);
}
