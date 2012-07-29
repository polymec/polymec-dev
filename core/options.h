#ifndef ARBI_OPTIONS_H
#define ARBI_OPTIONS_H

#include "arbi.h"

#ifdef __cplusplus
extern "C" {
#endif

// A type that stores command line options.
typedef struct options_t options_t;

// Parse options from the command line, returning an options object.
options_t* options_parse(int argc, char** argv);

// Free the storage allocated to the given options.
void options_free(options_t* opts);

// Returns the name of the computational model passed to the command line.
char* options_model(options_t* opts);

// Returns true if the command line requested help.
bool options_help(options_t* opts);

// Returns the name of any requested benchmark calculation, or NULL if no 
// such benchmark was requested.
char* options_benchmark(options_t* opts);

// Returns the input file passed to the command line, or NULL if no such 
// input file was given.
char* options_file(options_t* opts);

#ifdef __cplusplus
}
#endif

#endif

