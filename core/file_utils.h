// Copyright (c) 2012-2019, Jeffrey N. Johnson
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

/// \addtogroup files core:files
/// Functions for creating and manipulating files.
///@{

typedef struct string_slist_t string_slist_t;

/// Given a full pathname, parse it into directory and file portions.
/// Memory must be allocated for dirname and for filename that is sufficient
/// to store any portion of path.
/// \param [in] path An absolute path.
/// \param [out] dirname Stores the directory extracted from the path.
/// \param [out] filename Stores the filename extracted from the path.
void parse_path(const char* path, char* dirname, char* filename);

/// Given a path and a filename, join them using the OS-specific separator,
/// storing the result in path.
/// \param [in] dirname The absolute path to a directory.
/// \param [in] filename The name of a file (without any containing directory).
/// \param [out] path Stores the absolute path of the file within the directory.
void join_paths(const char* dirname, const char* filename, char* path);

/// Returns `true` if the given file exists (and is readable), `false` if not.
/// \param [in] filename The absolute path to a file.
bool file_exists(const char* filename);

/// Returns `true` if the given directory exists, `false` if not.
/// \param [in] dirname The absolute path to a directory.
bool directory_exists(const char* dirname);

/// Returns a linked list containing the names of all files within the
/// given directory. Filenames are relative to the directory. The directories
/// `.` and `..` are excluded from this list.
/// \param [in] dirname The absolute path to a directory.
string_slist_t* files_within_directory(const char* dirname);

/// Returns a linked list containing the names of all directories within
/// the given directory. Directory names are relative to the directory.
/// \param [in] dirname The absolute path to a directory.
string_slist_t* directories_within_directory(const char* dirname);

/// Creates the given directory if it doesn't exist and if it's capable.
/// If it does exist, this function has no effect.
/// \param [in] dirname An absolute path identifying a directory to create.
/// \param [in] mode The mode for the newly created directory. See mkdir (2)
///                  for details.
/// \returns `true` if the directory exists at the end of the call, `false`
/// otherwise.
bool create_directory(const char* dirname, mode_t mode);

/// Removes a directory and any of its contents.
/// \param [in] path An absolute path identifying a directory to remove.
/// \returns `true` on success, `false` on failure.
bool remove_directory(const char* path);

/// Creates a temporary file using the given template.
/// All temporary files and directories are deleted when a polymec application
/// exits. The temporary directory is created within a polymec-specific
/// directory unique to the running process.
/// \param [in] file_template A template containing 6 `X` characters to be replaced
///             with a set of characters that renders it unique (in the spirit
///             of `mkstemp`). Should look like this: `path/to/fileXXXXXX`, with the
///             `X`s all at the end.
/// \param [out] filename The absolute path to the created temporary file.
/// \returns A `FILE` object for the created file.
FILE* make_temp_file(const char* file_template, char* filename);

/// Creates a temporary directory using the given template.
/// All temporary files and directories are deleted when a polymec application
/// exits. The temporary directory is created within a polymec-specific
/// directory unique to the running process.
/// \param [in] dir_template A template containing 6 `X` characters to be replaced
///             with a set of characters that renders it unique (in the spirit
///             of `mkdtemp`). Should look like this: `path/to/dirXXXXXX`, with the
///             `X`s all at the end.
/// \param [out] dirname The absolute path to the created temporary directory.
/// \returns `true` if the directory was created, `false` if not.
bool make_temp_directory(const char* dir_template, char* dirname);

///@}
///@}

#endif
