/****
 * @date Created on 2021-04-06 at 12:54:13 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the real vector utilities.
 ***/
#ifndef _REAL_VECTOR_UTIL_H
#define _REAL_VECTOR_UTIL_H
#include <stdio.h>

void copy_vector(int n, double* src, double* dst);
double norm(int n, double* x);
double distance(int n, double* x, double* y);
void add_vector(int n, double* x, double* dx, double alpha);
void swap_vector(double* x, int i, int j);
void sort_vector(int n, double* x);
int uniq_vector(int* n, double* x, double toler);
double quantile(int n, double* x, double q);
double mean_vector(int n, double* x);
void print_vector(int n, double* x, const char* name);
void save_real_data(int n, double* x, const char* name, FILE* fp);
int count_real_data(const char* arg);
void parse_real_data(int n, double* data, char* arg);

#endif
