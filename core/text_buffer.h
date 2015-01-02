// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_TEXT_BUFFER_H
#define POLYMEC_TEXT_BUFFER_H

#include <stdbool.h>

// This is a line-by-line file buffer used for reading text files.
typedef struct text_buffer_t text_buffer_t;

// Creates a text buffer from the contents of the given file, or 
// returns NULL if the file cannot be opened.
text_buffer_t* text_buffer_from_file(const char* filename);

// Creates a text buffer from the given string.
text_buffer_t* text_buffer_from_string(const char* string);

// Destroys the file buffer.
void text_buffer_free(text_buffer_t* buffer);

// Returns the size (in bytes) of the buffer. 
long text_buffer_size(text_buffer_t* buffer);

// Returns the number of lines in the buffer.
int text_buffer_num_lines(text_buffer_t* buffer);

// Iterates over the file, returning the next line. Set pos to 0 
// to reset the iteration. line is not NULL-terminated--it is just the portion
// of the buffer holding the current line. Thus, you should 
// use line_length to determine its length.
bool text_buffer_next(text_buffer_t* buffer, int* pos, char** line, int* line_length);

// Iterates over the file, returning the next non-empty line. Set pos to 0 
// to reset the iteration. line is not NULL-terminated--it is just the portion
// of the buffer holding the current line. Thus, you should 
// use line_length to determine its length.
bool text_buffer_next_nonempty(text_buffer_t* buffer, int* pos, char** line, int* line_length);

// Returns a newly-allocated string containing the contents of this text buffer.
char* text_buffer_to_string(text_buffer_t* buffer);

#endif
