/****
 * @date Created on 2020-09-07 at 17:03:10 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the functions to parse configuration files.
 ***/
#ifndef _CONFIG_IO_H
#define _CONFIG_IO_H
#include "config_type.h"

Config* new_config(const char* fname);
void del_config(Config* conf);
void generate_filename(Config* conf, const char* outputdir, char* fname);

#endif
