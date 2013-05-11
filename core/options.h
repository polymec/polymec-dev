#ifndef POLYMEC_OPTIONS_H
#define POLYMEC_OPTIONS_H

#include "polymec.h"

// A type that stores command line options. Objects of this type are 
// garbage-collected.
typedef struct options_t options_t;

// Creates an empty options object. This is mostly for simplifying logic.
options_t* options_new();

// Parses options from the command line, returning an options object.
options_t* options_parse(int argc, char** argv);

// Returns the command passed on the command line, or NULL if no command 
// was given.
char* options_command(options_t* opts);

// Returns the input identifier passed to the command line, or NULL if no such 
// identifier was given. If the given command was "run", this identifies an
// input file. If the command was "benchmark", this identifies a problem.
char* options_input(options_t* opts);

// Returns a string value for the given parameter, or NULL if no such 
// parameter exists within opts.
char* options_value(options_t* opts, const char* name);

// Sets the given value for the the given option.
void options_set(options_t* opts, const char* name, const char* value);

#endif

