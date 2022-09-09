/****
 * @date Created on 2020-08-22 at 19:32:49 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C header providing the common purpose functions.
 ***/
#ifndef _COMMON_UTIL_H
#define _COMMON_UTIL_H

//int ipow(int x, int n); //Not used anywhere.
void replace_ext(const char* infname, const char* ext, char* outfname);
void split_path(const char* path, char* dir, char* base);
int is_char_in_str(const char c, const char* str);
int is_blank(const char c);
int is_space(const char c);
void strip_spaces(char* str);
int no_file(char* fname);
int no_text_file(const char* fname);
int ensure_path(const char* path);
void prepare_file(const char* fname, const char* comstr, const char* heading);
int convert_ppm_to_png(const char* fname);
int get_value_index(const char* str, const char* key);
char* get_value(int narg, char** args, const char* key);
void generate_uuid(char* uuid);
void print_progress(int inow, int ntot, double elapsed, const char* msg);
void print_end(int ntot, double elapsed, const char* msg);

#endif
