/****
 * @date Created on 2021-05-09 at 19:53:38 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing tag strings to highlight information with colors, in particular for tests.
 ***/
#ifndef _TAG_STRINGS_H
#define _TAG_STRINGS_H

#define STR_FAIL   "[\033[1;31mFAIL\033[0m] "  //String to highlight failed tests.
#define STR_PASS   "[\033[1;32mPASS\033[0m] "  //String to highlight successful tests.
#define STR_WARN   "[\033[1;33mWARN\033[0m] "  //String to highlight warnings.

#endif
