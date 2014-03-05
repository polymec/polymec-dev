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

#include <stdio.h>
#include "core/text_file_buffer.h"
#include "core/slist.h"

struct text_file_buffer_t 
{
  long size;
  int num_lines;
  long* line_offsets;
  char* data;
};

text_file_buffer_t* text_file_buffer_new(const char* filename)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL)
    return NULL;

  text_file_buffer_t* buffer = malloc(sizeof(text_file_buffer_t));

  // Read the contents of the entire file into memory.
  long start = ftell(fp);
  fseek(fp, 0, SEEK_END);
  long end = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  buffer->size = end - start;
  buffer->data = malloc(sizeof(char) * buffer->size);
  fread(buffer->data, sizeof(char), buffer->size,fp);

  // We're through with the file now.
  fclose(fp);

  // Now find all the line breaks. Here's why you don't want to use 
  // this for enormous files!
  long_slist_t* line_list = long_slist_new();
  for (long i = 0; i < buffer->size; ++i)
  {
    if (buffer->data[i] == '\n')
      long_slist_append(line_list, i+1);
  }

  buffer->num_lines = line_list->size;
  buffer->line_offsets = malloc(sizeof(long) * (buffer->num_lines+1));
  buffer->line_offsets[0] = 0;
  int i = 1;
  for (long_slist_node_t* iter = line_list->front; iter != NULL; iter = iter->next)
  {
    buffer->line_offsets[i] = iter->value;
    ++i;
  }
  long_slist_free(line_list);

  return buffer;
}

void text_file_buffer_free(text_file_buffer_t* buffer)
{
  free(buffer->data);
  free(buffer->line_offsets);
  free(buffer);
}

long text_file_buffer_size(text_file_buffer_t* buffer)
{
  return buffer->size;
}

int text_file_buffer_num_lines(text_file_buffer_t* buffer)
{
  return buffer->num_lines;
}

bool text_file_buffer_next(text_file_buffer_t* buffer, int* pos, char** line, int* line_length)
{
  if (*pos >= buffer->num_lines) 
    return false;
  int offset = buffer->line_offsets[*pos];
  *line_length = buffer->line_offsets[*pos + 1] - offset;
  *line = &buffer->data[offset];
  ++(*pos);
  return true;
}

