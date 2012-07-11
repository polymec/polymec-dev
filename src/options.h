#ifndef CLIDE_OPTIONS_H
#define CLIDE_OPTIONS_H

#include "clide.h"

#ifdef __cplusplus
extern "C" {
#endif

// A type that stores command line options.
typedef struct options_t options_t;

// Parse options from the command line.
options_t* options_parse(int argc, char** argv);

// Free the storage allocated to the given options.
void options_free(options_t* opts);

// Returns the input file passed to the command line, or NULL if no such 
// input file was given.
char* options_file(options_t* opts);

#ifdef __cplusplus
}
#endif

#endif

