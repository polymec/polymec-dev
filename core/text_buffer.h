// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_TEXT_BUFFER_H
#define POLYMEC_TEXT_BUFFER_H

#include <stdbool.h>

/// \addtogroup core core
///@{

/// \class text_buffer
/// This is a line-by-line file buffer used for reading text files.
typedef struct text_buffer_t text_buffer_t;

/// Creates a text buffer from the contents of the given file, or
/// returns NULL if the file cannot be opened.
/// \memberof text_buffer
text_buffer_t* text_buffer_from_file(const char* filename);

/// Creates a text buffer from the given string.
/// \memberof text_buffer
text_buffer_t* text_buffer_from_string(const char* string);

/// Destroys the file buffer.
/// \memberof text_buffer
void text_buffer_free(text_buffer_t* buffer);

/// Returns the size (in bytes) of the buffer.
/// \memberof text_buffer
size_t text_buffer_size(text_buffer_t* buffer);

/// Returns the number of lines in the buffer.
/// \memberof text_buffer
size_t text_buffer_num_lines(text_buffer_t* buffer);

/// Iterates over the file, returning the next line. Set pos to 0
/// to reset the iteration. line is not NULL-terminated--it is just the portion
/// of the buffer holding the current line. Thus, you should
/// use line_length to determine its length.
/// \memberof text_buffer
bool text_buffer_next(text_buffer_t* buffer, int* pos, char** line, size_t* line_length);

/// Iterates over the file, returning the next non-empty line. Set pos to 0
/// to reset the iteration. line is not NULL-terminated--it is just the portion
/// of the buffer holding the current line. Thus, you should
/// use line_length to determine its length.
/// \memberof text_buffer
bool text_buffer_next_nonempty(text_buffer_t* buffer, int* pos, char** line, size_t* line_length);

/// Returns a newly-allocated string containing the contents of this text buffer.
/// \memberof text_buffer
char* text_buffer_to_string(text_buffer_t* buffer);

///@}

#endif
