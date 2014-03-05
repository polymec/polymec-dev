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

#ifndef POLYMEC_TEXT_FILE_BUFFER_H
#define POLYMEC_TEXT_FILE_BUFFER_H

#include <stdbool.h>

// This is a line-by-line file buffer used for reading text files.
typedef struct text_file_buffer_t text_file_buffer_t;

// Creates a read-only file buffer examining the given file, or 
// returns NULL if the file cannot be opened.
text_file_buffer_t* text_file_buffer_new(const char* filename);

// Destroys the file buffer.
void text_file_buffer_free(text_file_buffer_t* buffer);

// Returns the size (in bytes) of the buffer. Yes, this is a 
// regular 32-bit integer. You shouldn't try to read files into memory
// whose sizes don't fit into an int, dummy.
int text_file_buffer_size(text_file_buffer_t* buffer);

// Returns the number of lines in the buffer.
int text_file_buffer_num_lines(text_file_buffer_t* buffer);

// Iterates over the file, returning the next line. Set pos to 0 
// to reset the iteration. Since line is not NULL-terminated, you should 
// use line_length to determine its length.
bool text_file_buffer_next(text_file_buffer_t* buffer, int* pos, char** line, int* line_length);

#endif
