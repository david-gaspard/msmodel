/***
 * @author David GASPARD <dgaspard@ulb.ac.be>
 * @date Created on 2020-08-03 at 19:16 CEST
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file Header file of the core source file providing the higher-levels functions.
 **/
#ifndef _MS_CORE_H
#define _MS_CORE_H
#include "config_io.h"

void ms_core_main(Config* conf, const char* outfname, int nthread);

#endif

