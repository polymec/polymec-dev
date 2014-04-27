// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#ifndef POLYMEC_FILE_UTILS_H
#define POLYMEC_FILE_UTILS_H

// Given a full pathname, parse it into directory and file portions.
// Memory must be allocated for dirname and for filename that is sufficient 
// to store any portion of path.
void parse_path(const char *path, char *dirname, char *filename);

// Given a path and a filename, join them using the OS-specific separator, 
// storing the result in path.
void join_paths(const char *dirname, const char* filename, char* path);

// Create a temporary file using the given template. This function replaces 
// all X characters in the filename template with a set of characters that 
// renders the filename unique (in the spirit of mkstemp). The template should
// have the form path/to/fileXXXXXX, with the X's all at the end. This 
// function returns a file descriptor that is open for writing ("w") on 
// success, or NULL on failure.
FILE* make_temp_file(char* filename);

// This version of make_temp_file creates a temporary file that is constructed 
// from the given template (filename), appending the given suffix to the end.
// It returns a file descriptor that is open for writing ("w") on success, or 
// NULL on failure.
FILE* make_temp_file_with_suffix(char* filename, const char* suffix);

#endif
