/****
 * @date Created on 2020-08-15 at 18:10:59 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the pixel type.
 ***/
#ifndef _PIXEL_TYPE_H
#define _PIXEL_TYPE_H
#include <stdint.h>

typedef struct Pixel_s {
	uint8_t red;
	uint8_t green;
	uint8_t blue;
} Pixel;

#endif
