/****
 * @date Created on 2020-08-22 at 19:17:54 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code for testing the general purpose functions.
 ***/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sys/stat.h>      /* Library providing the UNIX stat/mkdir Utility */
#include "common_util.h"   /* Import the General Purpose Functions */

/**
 * Test the replace_ext() function with some special cases.
 */
void test_replace_ext(void) {
	int i, fnamelen = 14;
	char fnamels[14][30] = {
		"myfile.txt",
		".emacs",
		"anotherfile.",
		"nodot",
		"",
		"dot.dot.dot",
		"/to/my/nodot",
		"/another/path.file/struct",
		"path/to/file.txt",
		"yet/another/.test",
		"double/dot/with..slash",
		"double..dot",
		"..ddot",
		"enddots..."
	};
	const char ext[] = ".png";
	char outfname[30];
	for (i = 0; i < fnamelen; i++) {
		printf("[TEST] File%d: %s", i, fnamels[i]);
		replace_ext(fnamels[i], ext, outfname);
		printf(" -> %s\n", outfname);
	}
}
/**
 * Test the split_path() function with some special cases.
 */
void test_split_path(void) {
	int i, pathlen = 11;
	char pathls[11][30] = {
		"myfile.txt",
		"",
		"path/to/my/file",
		"out/filename.txt",
		"~/path/to/f",
		"/usr",
		"/usr/",
		"/usr/lib",
		"~/usr/lib/",
		"/another/path.file/struct",
		"ill/formed/path//"
	};
	char dir[30], base[30];
	for (i = 0; i < pathlen; i++) {
		printf("[TEST] Path%d: %s", i, pathls[i]);
		split_path(pathls[i], dir, base);
		printf(" -> dir=%s, base=%s\n", dir, base);
	}
}

/**
 * Ensure that all directories in the given "path" exist. Creates the missing directories if needed.
 * The possible file at the end of "path" is omitted in the creation process, as long as it has no trailing slash.
 */
int mkpath(const char* path, mode_t mode) {
	//const char path[] = "d3/cube/a1.0/myfile.dat";
	const char sep = '/';
	char dir[strlen(path)+1];
	struct stat st;
	int created = 0; //Flag indicating whether a directory has been created.
	int i;
	for (i = 0; path[i]; i++) {
		if (path[i] == sep) {
			sprintf(dir, "%.*s", i+1, path);
			if (stat(dir, &st) != 0) {//If directory does not exist yet.
				if (mkdir(dir, mode) == 0) {//If directory properly created.
					created = 1;
				}
				else {
					printf("[ERROR] Failed to create directory '%s'.\n", dir);
					return 1;
				}
			}
		}
	}
	if (created) {
		printf("[INFO] Created directory for '%s'...\n", path);
	}
	return 0;
}

int main(int argc, char** argv) {
	
	//mkpath("yet/another/test/withmkpath/file.txt", 0700);
	//ensure_path("another/test/withensurepath/file.txt");
	
	//printf("[TEST] Testing extension replacement:\n");
	//test_replace_ext();
	
	//printf("[TEST] Testing path splitting:\n");
	//test_split_path();
	
	//char teststr1[] = "    To   be  or\t   not   to  be,       that\t\t\tis \t  the  question.   ";
	//printf("[INFO] Original string = '%s'\n", teststr1);
	//strip_spaces(teststr1);
	//printf("[INFO] Stripped string = '%s'\n", teststr1);
	//
	//char teststr2[] = "\nJust a normal phrase\tto   be sure.  ";
	//printf("[INFO] Original string = '%s'\n", teststr2);
	//strip_spaces(teststr2);
	//printf("[INFO] Stripped string = '%s'\n", teststr2);
	//
	//char teststr3[] = "";
	//printf("[INFO] Original string = '%s'\n", teststr3);
	//strip_spaces(teststr3);
	//printf("[INFO] Stripped string = '%s'\n", teststr3);
	//
	//char c = 'b';
	//if (is_char_in_str(c, teststr2)) {
	//	printf("[INFO] Character '%c' FOUND in '%s'\n", c, teststr2);
	//}
	//else {
	//	printf("[INFO] Character '%c' NOT found in '%s'\n", c, teststr2);
	//}
	
	//if (argc < 2) {
	//	printf("[USAGE] %s file\n", argv[0]);
	//	exit(1);
	//}
	//char* fname = argv[1];
	//if (no_file(fname)) {
	//	printf("[WARN] No file '%s' found...\n", fname);
	//}
	//else {
	//	printf("[INFO] OK, file '%s' exists. We can proceed...\n", fname);
	//}
	//if (no_text_file(fname)) {
	//	printf("[WARN] The given file '%s' is NOT a text file...\n", fname);
	//}
	//else {
	//	printf("[INFO] OK, file '%s' is a text file...\n", fname);
	//}
	
	return 0;
}
