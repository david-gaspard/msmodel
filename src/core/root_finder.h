/****
 * @date Created on 2021-04-18 at 18:02:46 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing some root-finding algorithms for complex functions.
 ***/
#ifndef _ROOT_FINDER_H
#define _ROOT_FINDER_H
#include "domain_type.h"

dcomplex eval_cpoly(int d, dcomplex* poly, dcomplex z);
void fit_poly(int ns, double* xs, double* fs, int d, double* poly);
void fit_cpade(int ns, dcomplex* zs, dcomplex* fs, int d, dcomplex* pade);
void find_root_poly(int d, dcomplex* poly, dcomplex* roots);
int find_root_secant(dcomplex (*f)(dcomplex), Domain* dom, dcomplex* z, int maxit, double toler, int verb);
int find_root_muller(dcomplex (*f)(dcomplex), Domain* dom, dcomplex* z, int maxit, double toler, int verb);
void guess_roots(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int verb);
void solve_from_guess(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int maxit, double toler, int verb);
void solve_aberth(dcomplex (*dlogf)(dcomplex), Domain* dom, int* nz, dcomplex* z, int nexpsup, int maxit, double toler, int verb);
void solve_aberth_num(dcomplex (*f)(dcomplex), Domain* dom, int* nz, dcomplex* z, int maxit, double toler, int verb);
int find_fit_gna(double (*f)(int i, double* p), int ns, double* ys, int np, double* param, int maxit, int verb);

#endif
