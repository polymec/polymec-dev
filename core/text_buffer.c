// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <stdio.h>
#include "core/text_buffer.h"
#include "core/slist.h"

struct text_buffer_t 
{
  long size;
  int num_lines;
  long* line_offsets;
  char* data;
};

static void find_line_breaks(text_buffer_t* buffer)
{
  // Here's why you don't want to use this for enormous files!
  long_slist_t* line_list = long_slist_new();
  for (long i = 0; i < buffer->size; ++i)
  {
    if (buffer->data[i] == '\n')
      long_slist_append(line_list, i+1);
  }
  {
    // Pick up all the detritus at the end if necessary.
    long_slist_node_t* iter = line_list->front;
    long last;
    while (long_slist_next(line_list, &iter, &last));
    if (last < buffer->size)
      long_slist_append(line_list, buffer->size);
  }

  buffer->num_lines = line_list->size;
  buffer->line_offsets = polymec_malloc(sizeof(long) * (buffer->num_lines+1));
  buffer->line_offsets[0] = 0;
  int i = 1;
  long line_no;
  long_slist_node_t* iter = NULL;
  while (long_slist_next(line_list, &iter, &line_no))
  {
    buffer->line_offsets[i] = line_no; //iter->value;
    ++i;
  }
  long_slist_free(line_list);
}

text_buffer_t* text_buffer_from_file(const char* filename)
{
  FILE* fp = fopen(filename, "r");
  if (fp == NULL)
    return NULL;

  text_buffer_t* buffer = polymec_malloc(sizeof(text_buffer_t));

  // Read the contents of the entire file into memory.
  long start = ftell(fp);
  fseek(fp, 0, SEEK_END);
  long end = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  buffer->size = end - start + 1;
  buffer->data = polymec_malloc(sizeof(char) * buffer->size);
  size_t bytes_read = fread(buffer->data, sizeof(char), buffer->size, fp);
  buffer->data[bytes_read-1] = '\0';

  // We're through with the file now.
  fclose(fp);

  // Now find all the line breaks.
  find_line_breaks(buffer);

  return buffer;
}

text_buffer_t* text_buffer_from_string(const char* string)
{
  text_buffer_t* buffer = polymec_malloc(sizeof(text_buffer_t));

  // Copy the string.
  buffer->size = strlen(string) + 1;
  buffer->data = polymec_malloc(sizeof(char) * buffer->size);
  memcpy(buffer->data, string, sizeof(char) * buffer->size);

  // Now find all the line breaks.
  find_line_breaks(buffer);

  return buffer;
}

void text_buffer_free(text_buffer_t* buffer)
{
  polymec_free(buffer->data);
  polymec_free(buffer->line_offsets);
  polymec_free(buffer);
}

long text_buffer_size(text_buffer_t* buffer)
{
  return buffer->size;
}

int text_buffer_num_lines(text_buffer_t* buffer)
{
  return buffer->num_lines;
}

bool text_buffer_next(text_buffer_t* buffer, int* pos, char** line, int* line_length)
{
  if (*pos >= buffer->num_lines) 
    return false;
  int offset = (int)buffer->line_offsets[*pos];
  *line_length = (int)buffer->line_offsets[*pos + 1] - offset - 1;
  *line = &buffer->data[offset];
  ++(*pos);
  return true;
}

bool text_buffer_next_nonempty(text_buffer_t* buffer, int* pos, char** line, int* line_length)
{
  bool result = text_buffer_next(buffer, pos, line, line_length);
  while (result && (line_length == 0))
    result = text_buffer_next(buffer, pos, line, line_length);
  return result;
}

char* text_buffer_to_string(text_buffer_t* buffer)
{
  return string_dup(buffer->data);
}

