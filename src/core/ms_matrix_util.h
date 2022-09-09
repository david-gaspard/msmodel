/****
 * @date Created on 2020-08-09 at 14:03:47 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the tools to create and manage the multiple scattering (MS) matrix.
 ***/
#ifndef _MS_MATRIX_UTIL_H
#define _MS_MATRIX_UTIL_H
#include "dcomplex_type.h"
#include "medium_type.h"

void build_ms_matrix(Medium* med, dcomplex k, dcomplex* matrix);
void build_normalized_matrix(Medium* med, dcomplex k, dcomplex* matrix);
void build_dk_ms_matrix(Medium* med, dcomplex k, dcomplex* matrix);
void build_plane_wave(Medium* med, dcomplex k, double theta, dcomplex* vec);
void print_matrix(int n, dcomplex* matrix, int s);
void eigvals(int n, dcomplex* matrix, dcomplex* mu);
double tracegtwo(int n, dcomplex* matrix);
dcomplex tracelog_lu(int n, dcomplex* matrix);
//dcomplex tracelog_ev(int n, dcomplex* matrix);  //Deprecated for inefficiency.
//dcomplex tracelog_qr(int n, dcomplex* matrix);  //Deprecated for inefficiency (not more stable than LU).
dcomplex dk_log_det_ms(Medium* med, dcomplex k);
int inverse_iteration(int n, dcomplex* a, dcomplex* vec, dcomplex* lambda, int maxit, double toler, int verb);
dcomplex traceinv(int n, dcomplex* matrix);
double total_cross_section(Medium* med, double k);

#endif
