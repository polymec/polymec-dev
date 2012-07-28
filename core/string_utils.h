#ifndef ARBI_STRING_UTILS_H
#define ARBI_STRING_UTILS_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// Given a full pathname, parse it into directory and file portions.
// Memory must be allocated for dirname and for filename that is sufficient 
// to store any portion of path.
int parse_path(const char *path, char *dirname, char *filename);

// Given a path and a filename, join them using the OS-specific separator, 
// storing the result in path.
int join_paths(const char *dirname, const char* filename, char* path);

#ifdef __cplusplus
}
#endif

#endif
