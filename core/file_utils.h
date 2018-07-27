// Copyright (c) 2012-2018, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_FILE_UTILS_H
#define POLYMEC_FILE_UTILS_H

#include <sys/types.h> // for mode_t.
#include <sys/stat.h> // for mkdir().

/// \addtogroup core core
///@{

typedef struct string_slist_t string_slist_t;

/// Given a full pathname, parse it into directory and file portions.
/// Memory must be allocated for dirname and for filename that is sufficient 
/// to store any portion of path.
void parse_path(const char* path, char* dirname, char* filename);

/// Given a path and a filename, join them using the OS-specific separator, 
/// storing the result in path.
void join_paths(const char* dirname, const char* filename, char* path);

/// Returns true if the given file exists (and is readable), false if not.
bool file_exists(const char* filename);

/// Returns true if the given directory exists, false if not.
bool directory_exists(const char* dirname);

/// Returns a linked list containing the names of all files within the 
/// given directory. Filenames are relative to the directory. The directories
/// . and .. are excluded from this list.
string_slist_t* files_within_directory(const char* dirname);

/// Returns a linked list containing the names of all directories within 
/// the given directory. Directory names are relative to the directory.
string_slist_t* directories_within_directory(const char* dirname);

/// Creates the given directory if it doesn't exist and if it's capable. 
/// If it does exist, this function has no effect. Returns true if the 
/// directory exists at the end of the call, false otherwise.
bool create_directory(const char* dirname, mode_t mode);

/// Remove a directory and any of its contents. Returns 0 on success, -1 on 
/// failure (as does the standard rmdir).
int remove_directory(const char* path);

/// Create a temporary file using the given template. This function replaces 
/// 6 X characters in the filename template with a set of characters that 
/// renders it unique (in the spirit of mkstemp), storing the result in filename. 
/// The template should have the form path/to/fileXXXXXX, with the X's all at 
/// the end. This function returns a file descriptor that is open for writing 
/// ("w") on success, or NULL on failure. All temporary files are deleted 
/// when a polymec application exits. Note that the temporary file 
/// is created within a unique polymec-specific temporary directory.
FILE* make_temp_file(const char* file_template, char* filename);

/// Create a temporary directory using the given template. This function replaces 
/// 6 X characters in the dirname template with a set of characters that 
/// renders it unique (in the spirit of mkdtemp). The template should have the 
/// form path/to/dirXXXXXX, with the X's all at the end. This function returns 
/// true if the directory was created, false if not. All temporary files are 
/// deleted when a polymec application exits. Note that the temporary directory 
/// is created within a unique polymec-specific temporary directory.
bool make_temp_directory(const char* dir_template, char* dirname);

///@}

#endif
