// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef POLYMEC_STRING_UTILS_H
#define POLYMEC_STRING_UTILS_H

// Given a full pathname, parse it into directory and file portions.
// Memory must be allocated for dirname and for filename that is sufficient 
// to store any portion of path.
void parse_path(const char *path, char *dirname, char *filename);

// Given a path and a filename, join them using the OS-specific separator, 
// storing the result in path.
void join_paths(const char *dirname, const char* filename, char* path);

// Since strdup() is not standard C, we provide a surrogate here.
char* string_dup(const char* s);

// Given a string and a NULL-terminated list of token-value pairs, 
// returns a newly-allocated string containing the original string with 
// all instances of tokens replaced with their substitution values. This 
// string must be freed with free().
typedef struct { char* token; char* value; } string_subst_t;
static const string_subst_t END_OF_SUBST = {(char*)"END_TOKEN", (char*)"END_VALUE"};
char* string_subst(const char* string, string_subst_t substitutions[]);

#endif
