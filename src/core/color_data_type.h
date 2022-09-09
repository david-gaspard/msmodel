/****
 * @date Created on 2021-03-21 at 12:09:02 CET
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the color data type which provides a pixel for any real between 0 and 1.
 ***/
#ifndef _COLOR_DATA_TYPE_H
#define _COLOR_DATA_TYPE_H

typedef struct ColorData_s {
	char name[50];              //Name of the color rule.
	int len;                    //Number of samples in the color data.
	const double (*colors)[3];  //Pointer to the first element of the color data.
} ColorData;

#endif
