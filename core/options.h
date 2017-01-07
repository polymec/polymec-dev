// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_OPTIONS_H
#define POLYMEC_OPTIONS_H

#include "core/polymec.h"

// A type that stores command line options. Objects of this type are 
// garbage-collected.
typedef struct options_t options_t;

// Creates an empty options object. This is mostly for simplifying logic.
options_t* options_new(void);

// Returns the options_argv singleton that holds all of the parsed command 
// line options.
options_t* options_argv(void);

// Parses options from the command line, initializing the command line 
// options singleton.
void options_parse(int argc, char** argv);

// Returns an internal string holding the numbered argument passed on the 
// command line, in the same way argv would store it, or NULL if no such 
// argument was given.
char* options_argument(options_t* opts, int n);

// Returns the number of arguments given on the command line.
int options_num_arguments(options_t* opts);

// Returns a string holding the given named value, or NULL if no such 
// parameter exists within opts.
char* options_value(options_t* opts, const char* name);

// Sets the given value for the the given option.
void options_set(options_t* opts, const char* name, const char* value);

#endif

