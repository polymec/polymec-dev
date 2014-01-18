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

#ifndef POLYMEC_STRING_UTILS_H
#define POLYMEC_STRING_UTILS_H

#include <stdbool.h>

// Given a full pathname, parse it into directory and file portions.
// Memory must be allocated for dirname and for filename that is sufficient 
// to store any portion of path.
void parse_path(const char *path, char *dirname, char *filename);

// Given a path and a filename, join them using the OS-specific separator, 
// storing the result in path.
void join_paths(const char *dirname, const char* filename, char* path);

// Since strdup() is not standard C, we provide a surrogate here.
char* string_dup(const char* s);

// Returns the number of substrings (separated by the given delimiter)
// occur in the given string.
int string_num_tokens(const char* s, const char* delimiter);

// Split a string into substrings using the given delimiter, and return 
// an array of newly-allocated strings. num_substrings will contain the 
// length of this array. If the delimitor is not found, the whole string will 
// be a single substring.
char** string_split(const char* s, const char* delimiter, int* num_substrings);

// Returns true if the given string is numeric, false if not.
bool string_is_number(const char* s);

// Given a string and a NULL-terminated list of token-value pairs, 
// returns a newly-allocated string containing the original string with 
// all instances of tokens replaced with their substitution values. This 
// string must be freed with free().
typedef struct { char* token; char* value; } string_subst_t;
static const string_subst_t END_OF_SUBST = {(char*)"END_TOKEN", (char*)"END_VALUE"};
char* string_subst(const char* string, string_subst_t substitutions[]);

#endif
