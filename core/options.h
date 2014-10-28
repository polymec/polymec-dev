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

#ifndef POLYMEC_OPTIONS_H
#define POLYMEC_OPTIONS_H

#include "core/polymec.h"

// A type that stores command line options. Objects of this type are 
// garbage-collected.
typedef struct options_t options_t;

// Creates an empty options object. This is mostly for simplifying logic.
options_t* options_new();

// Returns the options_argv singleton that holds all of the parsed command 
// line options.
options_t* options_argv();

// Parses options from the command line, initializing the command line 
// options singleton.
void options_parse(int argc, char** argv);

// Returns an internal string holding the numbered argument passed on the 
// command line, in the same way argv would store it, or NULL if no such 
// argument was given.
char* options_argument(options_t* opts, int n);

// Returns a string holding the given named value, or NULL if no such 
// parameter exists within opts.
char* options_value(options_t* opts, const char* name);

// Sets the given value for the the given option.
void options_set(options_t* opts, const char* name, const char* value);

#endif

