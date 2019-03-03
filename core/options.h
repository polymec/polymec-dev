// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_OPTIONS_H
#define POLYMEC_OPTIONS_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// \class options
/// A type that stores command line options.
/// \refcounted
typedef struct options_t options_t;

/// Creates an empty options object. This is mostly for simplifying logic.
/// \memberof options
options_t* options_new(void);

/// Returns the options_argv singleton that holds all of the parsed command
/// line options.
/// \memberof options
options_t* options_argv(void);

/// Parses options from the command line, initializing the command line
/// options singleton.
/// \memberof options
void options_parse(int argc, char** argv);

/// Returns an internal string holding the numbered argument passed on the
/// command line, in the same way argv would store it, or NULL if no such
/// argument was given.
/// \memberof options
char* options_argument(options_t* opts, size_t n);

/// Returns the number of arguments given on the command line, including
/// parameters with values.
/// \memberof options
size_t options_num_arguments(options_t* opts);

/// Returns true if the given argument was passed to the command line,
/// with or without an associated value. False otherwise.
/// \memberof options
bool options_has_argument(options_t* opts, const char* arg);

/// Appends the given argument to the list.
/// \memberof options
void options_add_argument(options_t* opts, const char* arg);

/// Removes the given argument from the list. Use with caution. This has
/// no effect if the given argument doesn't exist.
/// \memberof options
void options_remove_argument(options_t* opts, size_t n);

/// Returns a string holding the given named value, or NULL if no such
/// parameter exists within opts.
/// \memberof options
char* options_value(options_t* opts, const char* name);

/// Sets the given value for the the given option.
/// \memberof options
void options_set(options_t* opts, const char* name, const char* value);

/// Use this to traverse the named values in the list. Set *pos to 0 to reset.
/// \memberof options
bool options_next_value(options_t* opts, int* pos, const char** name, const char** value);

///@}

#endif

