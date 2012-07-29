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

// Returns the command passed on the command line, or NULL if no command 
// was given.
char* options_command(options_t* opts);

// Returns the input identifier passed to the command line, or NULL if no such 
// identifier was given. If the given command was "run", this identifies an
// input file. If the command was "benchmark", this identifies a problem.
char* options_input(options_t* opts);

#ifdef __cplusplus
}
#endif

#endif

