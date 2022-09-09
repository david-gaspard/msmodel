/****
 * @date Created on 2021-03-12 at 16:47:02 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the complex domain type.
 ***/
#ifndef _DOMAIN_TYPE_H
#define _DOMAIN_TYPE_H
#include "dcomplex_type.h"  //Since this is a complex domain, the complex library will necessarily be used.

/**
 * Defines a domain in the complex plane.
 */
typedef struct Domain_s {
	double xmin;  //Minimum real part of "z".
	double xmax;  //Maximum real part of "z".
	double ymin;  //Minimum imaginary part of "z".
	double ymax;  //Maximum imaginary part of "z".
} Domain;

#endif
