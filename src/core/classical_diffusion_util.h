/****
 * @date Created on 2021-02-26 at 11:14:59 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the utilities to solve the classical diffusion problem.
 ***/
#ifndef _CLASSICAL_DIFFUSION_UTIL_H
#define _CLASSICAL_DIFFUSION_UTIL_H

double helmholtz_series(int d, double k, double r);
double dr_helmholtz(int d, double k, double r);
double dk_helmholtz(int d, double k, double r);
double dkdr_helmholtz(int d, double k, double r);
int diffusion_result_plot_v1(int d, double radius, double nsigma, char* str);
int diffusion_result_plot_v2(int d, double radius, double nsigma, char* str);

#endif
