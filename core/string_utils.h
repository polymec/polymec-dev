#ifndef POLYMEC_STRING_UTILS_H
#define POLYMEC_STRING_UTILS_H

#include "polymec.h"

// Given a full pathname, parse it into directory and file portions.
// Memory must be allocated for dirname and for filename that is sufficient 
// to store any portion of path.
void parse_path(const char *path, char *dirname, char *filename);

// Given a path and a filename, join them using the OS-specific separator, 
// storing the result in path.
void join_paths(const char *dirname, const char* filename, char* path);

// Given a string and a NULL-terminated list of token-value pairs, 
// returns a newly-allocated string containing the original string with 
// all instances of tokens replaced with their substitution values. This 
// string must be freed with free().
typedef struct { char* token; char* value; } string_subst_t;
static const string_subst_t END_OF_SUBST = {"END_TOKEN", "END_VALUE"};
char* string_subst(const char* string, string_subst_t substitutions[]);

#endif
