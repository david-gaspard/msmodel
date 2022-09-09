/****
 * @date Created on 2020-08-22 at 17:58:58 CEST
 * @author David Gaspard <dgaspard@ulb.ac.be>
 * @copyright Copyright (C) 2022  David Gaspard
 * @license This program is free software; it can be redistributed and/or modified under
 * the terms of the GNU General Public License v3.0 as published by the Free Software Foundation.
 * This program is distributed in the hope that it will be useful but WITHOUT ANY WARRANTY
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * The license file, titled "LICENSE", can be found in the root directory of this project.
 * @file C code providing common purpose functions such as filename manipulation, UUID generation, progress bar, etc.
 ***/
#include <stdlib.h>     /* Provides the exit function*/
#include <stdio.h>      /* Standard Library of Input and Output */
#include <string.h>     /* Standard Library of String Manipulation */
#include <time.h>       /* Standard Library for Time and Date */
#include <math.h>       /* Standard Library of Mathematical Functions */
#include <unistd.h>     /* Import the UNIX Standard Library to get Host/User Name */
#include <sys/stat.h>   /* Library providing the UNIX stat/mkdir Utility */

/***********************************
 * GENERAL MATHEMATICAL FUNCTIONS:
 **********************************/
//See: green_bessel.c
/**
 * Fast integer power for positive powers without error/overflow checking.
 */
int ipow(int x, int n) {
	int res = 1;
	while (n != 0) {
		if (n & 1) {
			res *= x;
		}
		x *= x;
		n >>= 1;
	}
	return res;
}

/********************
 * REAL VECTOR TOOLS
 *******************/
//See: real_vector_util.c

/***********************
 * COMPLEX VECTOR TOOLS
 **********************/
//See: complex_vector_util.c

/********************************
 * BASIC FILENAME MANIPULATIONS:
 *******************************/
/**
 * Replace the extension in the given filename "fname" by the given one "ext" and outputs the result to "outfname", assuming Posix file path convention.
 * The string "infname" is not altered in the process. The extension should contain the final dot.
 * An empty extension "ext" will result in a extension-stripped filename "outfname".
 * The output string should be initialized with char outfname[strlen(infname)+strlen(ext)+1] to ensure the correct size.
 */
void replace_ext(const char* infname, const char* ext, char* outfname) {
	int len = strlen(infname)-1;
    while (len > 0 && infname[len] != '.' && infname[len] != '/') {
        len--;
    }
    if (len == 0 || infname[len] == '/' || infname[len-1] == '/') {
		len = strlen(infname);
	}
    sprintf(outfname, "%.*s%s", len, infname, ext);
}
/**
 * Split the given file path "path" into the directory part "dir" and the base name "base", assuming Posix file path convention.
 * By definition, the last character of "dir" is '/', and the filename "base" does never begin with '/',
 * so that the direct concatenation "dir"+"base" is guaranteed to be meaningful and conforms to the original "path".
 * Note that "base" can also be a directory if it ends with '/'.
 */
void split_path(const char* path, char* dir, char* base) {
	int i = strlen(path)-1;
	if (i < 0) {//Protects against zero length path.
		i = 0;
	}
	while (i > 0 && path[i-1] != '/') {
		i--;
	}
    if (i == 0) {//If no directory part found, then assume the current directory.
		strcpy(dir, "./");
	}
	else {
		strncpy(dir, path, i);
		dir[i] = '\0';
	}
    strcpy(base, path+i);
}
/**
 * Searches for the presence of the character "c" in the array "str", omitting the trailing null character.
 * Returns 1 if "c" is found, or 0 if it is not.
 */
int is_char_in_str(const char c, const char* str) {
	int i;
	for (i = 0; str[i]; i++) {
		if (c == str[i]) {
			return 1;
		}
	}
	return 0;
}
/**
 * Detect white space characters, but NOT line return.
 */
int is_blank(const char c) {
	return is_char_in_str(c, " \t\v"); // c == ' ' || c == '\t' || c == '\v';
}
/**
 * Returns 1 if the character "c" is a space-like character (space, tab, line return, line feed), or zero otherwise.
 */
int is_space(const char c) {
	return is_char_in_str(c, " \t\v\n\f"); // c == ' ' || c == '\n' || c == '\t' || c == '\v' || c == '\f';
}
/**
 * Strip extra space-like characters from the given string.
 * Adapted from: https://stackoverflow.com/questions/17770202/remove-extra-whitespace-from-a-string-in-c
 */
void strip_spaces(char* str) {
	int i, x = 0;
	for (i = 0; str[i]; i++) {
		if (!is_space(str[i]) || (i > 0 && !is_space(str[i-1]))) {
			if (is_space(str[i])) {//Convert space-like character to regular space.
				str[x] = ' ';
			}
			else {
				str[x] = str[i];
			}
			x++;
		}
	}
	if (x > 0 && is_space(str[x-1])) {//Remove trailing space.
		x--;
	}
	str[x] = '\0';
}
/**
 * Check that a file "fname" exists and is readable.
 * Returns 0 if the file "fname" exists, or 1 if it does not (error flag).
 */
int no_file(char* fname) {
	return access(fname, F_OK) != 0;
}
/**
 * Check if the given file "fname" is an ASCII text file.
 * Returns 0 if the file "fname" is a text file, or 1 if it is not (error flag).
 * Note that 1 is also returned if the file does not exist, so there is no need for two checks.
 */
int no_text_file(const char* fname) {
	if (system("which file >/dev/null")) {
		printf("[WARN] Command 'file' does not exist.\n");
		return 1; //Detect an issue by default for safety.
	}
	char cmd[50+strlen(fname)];
	sprintf(cmd, "file %s | grep -wi \"text\\|ascii\" >/dev/null", fname);
	return system(cmd) != 0;
}
/**
 * Ensure that all directories in the given "path" exist. Creates the missing directories with permissions 0700 (ugo=rwx|---|---) if needed.
 * The possible file at the end of "path" is omitted in the creation process, as long as it has no trailing slash.
 */
int ensure_path(const char* path) {       
	mode_t old_mask = umask(0);  //Removes the previous file mode mask.
	const mode_t mode = 0700;    //Most private directory permissions (rwx|---|---).
	const char sep = '/';        //UNIX separator.
	char dir[strlen(path)+1];
	struct stat st;
	int i;
	for (i = 0; path[i]; i++) {
		if (path[i] == sep) {
			sprintf(dir, "%.*s", i+1, path);
			if (stat(dir, &st) != 0) {//If directory does not exist yet.
				if (mkdir(dir, mode) != 0) {//If directory no properly created.
					printf("[ERROR] Failed to create directory '%s'.\n", dir);
					umask(old_mask); //Reset the previous mask.
					return 1;
				}
			}
		}
	}
	umask(old_mask); //Reset the previous mask.
	return 0;
}
/**
 * Create a new file "fname" with permission 0600 (ugo=rw-|---|---), and insert a short header timestamp at the beginning.
 * The "heading" should contain the comment character, such as "#" or "%", followed by a past passive participle, such as "Created".
 * Typical format: #Created 2020-08-22 18:08:32 CEST in MSModel by <username@cluster.ulb.ac.be>
 */
void prepare_file(const char* fname, const char* comstr, const char* heading) {
	const mode_t mode = 0600; //Most private file permission for files (rw-|---|---).
	const size_t buflen = 50;
	time_t t;
    time(&t);
	struct tm* timeinfo = localtime(&t);
	char datestr[buflen], username[buflen], hostname[buflen];
	strftime(datestr, buflen, "%F %T %Z", timeinfo);
	if (gethostname(hostname, buflen) != 0) {//If hostname not found
		strcpy(hostname, "???");
	}
	if (getlogin_r(username, buflen) != 0) {//If username not found
		strcpy(username, "???");
	}
	FILE* fp = fopen(fname, "w");
	fprintf(fp, "%s%s %s in MSModel by <%s@%s>\n%sOriginal filename: %s\n", comstr, heading, datestr, username, hostname, comstr, fname);
	fclose(fp);
	chmod(fname, mode);
	//printf("[INFO] Changed file permission of '%s' to %04o\n", fname, mode);
}

/**
 * Convert a Portable PixMap (PPM) file "fname.ppm" to a PNG file "fname.png".
 * The filename "fname" should not have any extension.
 * Returns 0 on success, and 1 on failure.
 */
int convert_ppm_to_png(const char* fname) {
	char ppmfile[strlen(fname)+5];
	sprintf(ppmfile, "%s.ppm", fname);
	if (no_file(ppmfile)) {
		printf("[ERROR] File '%s' not found.\n", ppmfile);
		return 1;
	}
	const char converter[] = "pnmtopng";
	char test[30+strlen(converter)];
	sprintf(test, "which %s > /dev/null 2>&1", converter);
	if (system(test)) {
		printf("[ERROR] Command '%s' does not exist. No conversion PPM -> PNG was made.\n", converter);
		return 1;
	}
	char cmd[50+3*strlen(fname)];
	sprintf(cmd, "pnmtopng %s > %s.png 2>/dev/null; rm %s", ppmfile, fname, ppmfile);
	system(cmd);
	return 0;
}

/*****************************************************
 * PARSING FUNCTIONS AND KEY/VALUE ARRAY MANIPULATION
 ****************************************************/
/**
 * Detect if the given string "str" begins with the "key". Take account of the key/value separator.
 * Returns the index of the value following the "key" in "str".
 * For instance, str="dimension3" and key="dimension" should not match, but str="dimension=3" will.
 */
int get_value_index(const char* str, const char* key) {
	int i;
	for (i = 0; key[i]; i++) {
		if (str[i] != key[i]) {
			return 0;
		}
	}
	if (str[i] == '=' || str[i] == ':') {//Check for separator.
		return i+1;
	}
	return 0; //No separator.
}

/**
 * Finds the value in the argument list "args" corresponding to the given "key".
 */
char* get_value(int narg, char** args, const char* key) {
	int i, l;
	for (l = 0; l < narg; l++) {//Loop on arguments.
		i = get_value_index(args[l], key);
		if (i != 0) {//Compare the first characters of "key" and "args[l]" (including separator).
			return args[l]+i;
		}
	}
	return ""; //If key not found, then returning empty string.
}

/********************
 * UUID GENERATION:
 *******************/
/**
 * List of base characters in order:
 */
static const char BASESTR[] = "0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz@#"; //Length: 10+2*26+2 = 64.
/**
 * Converts the given number "val" (long int, 64-bit) to a string of given base. Writes the resulting string to "res".
 * The length of the car array of "res" can be smaller or equal to 64.
 */
void inttostr(char* res, long int val, size_t base) {
	size_t maxbase = strlen(BASESTR);
	if (base > maxbase) {
		printf("[WARN] The base is too large and has been set to the maximum base %lu\n", maxbase);
		base = maxbase;
	}
	int i, len = 0;
	while (val > 0) {//Base conversion (modular method)
		res[len] = BASESTR[val%base];
		val /= base;
		len++;
	}
	char tmp;
	for (i = 0; i < len/2; i++) {//Reverse the string.
		tmp = res[i];
		res[i] = res[len-i-1];
		res[len-i-1] = tmp;
	}
	res[len] = '\0'; //Null character at the end.
}
/**
 * Exploits the current time to generate a unique ID for the simulations.
 * More precisely, the number of milliseconds since the UNIX Epoch in UTC time zone (1970-01-01 at 00:00:00 UTC)
 * is converted in base 36 using the uppercase letters of the latin alphabet as the symbols beyond '9'.
 * The resulting ID is written to the char array "uuid". This array can be about 20 length.
 * The UUID can be used in the format 'file_uuid.tab' to name files. Sample UUID: 'KDSTQ67H'
 */
void generate_uuid(char* uuid) {
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);  //Get UNIX timestamp up to nanosecond resolution.
	if (now.tv_sec < 1600000000) {//Check if the UNIX Epoch Time is plausible given we are beyond 2020-09-14 right now.
		printf("[WARN] Inconsistent UNIX timestamp %ld s. The UUID will be offset.\n", now.tv_sec);
	}
	inttostr(uuid, now.tv_sec*1000 + now.tv_nsec/1000000, 36); //Convert the UNIX Epoch Time (in milliseconds) to base 36.
	//printf("[INFO] Generated UUID '%s'.\n", uuid);
}

/*************************
 * PROGRESS BAR DISPLAY:
 ************************/
/**
 * If zero number of days, then prints the given time in seconds in the standard format hh:mm:ss (without line break).
 * If nonzero number of days, then prints the given time in the format 0d 0h 0m 0s (without line break).
 */
void print_time(double time_in_sec) {
	int days = (int)floor(time_in_sec/86400);
	int hours = (int)floor(fmod(time_in_sec, 86400)/3600);
	int minutes = (int)floor(fmod(time_in_sec, 3600)/60);
	int seconds = (int)floor(fmod(time_in_sec, 60));
	if (days == 0) {//Standard time format (ISO-8601):
		printf("%02d:%02d:%02d", hours, minutes, seconds);
	}
	else {//More verbose format:
		printf("%2dd %02dh %02dm %02ds", days, hours, minutes, seconds);
	}
}
/**
 * Print an ASCII progress bar with the current status and the time left. This function erases the last line of the console.
 * "inow" is the current number of jobs done (should be nonzero), "ntot" is the total number of jobs, and
 * "elapsed" is the elapsed time since the beginning (inow=0) in seconds.
 * FORMAT: [EXEC] <msg_here> | 00:00:00 [###################______] 99.99% ETA 00:00:00
 */
void print_progress(int inow, int ntot, double elapsed, const char* msg) {
	const int len = 30; //Full length of the progress bar in characters (improvement: detect console width).
	double time_left = (ntot - inow)*elapsed/inow; //Estimated Time of Arrival (ETA).
	double progress = (double)inow/ntot;
	int pos = (int)round(progress*len);  //Current position of the progress bar.
	int i;
	printf("\r[EXEC] %s | ", msg);
	print_time(elapsed);
	printf(" [");
	for (i = 0; i < len; i++) {//Print the progress bar.
		if (i < pos) {
			printf("#");
		}
		else {
			printf("_");
		}
	}
	printf("] %5.2f%% ETA ", 100.0*progress);
	print_time(time_left);
	fflush(stdout);
}
/**
 * Prints the ending line of the progress bar with the total amount of time in seconds.
 */
void print_end(int ntot, double elapsed, const char* msg) {
	print_progress(ntot, ntot, elapsed, msg);
	printf("... Done\n");
	//printf("\n[DONE] All jobs completed.\n");
}
