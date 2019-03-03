// Copyright (c) 2012-2019, Jeffrey N. Johnson
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
  size_t size;
  size_t num_lines;
  size_t* line_offsets;
  char* data;
};

static void find_line_breaks(text_buffer_t* buffer)
{
  // Here's why you don't want to use this for enormous files!
  size_t_slist_t* line_list = size_t_slist_new();
  for (size_t i = 0; i < buffer->size; ++i)
  {
    if (buffer->data[i] == '\n')
      size_t_slist_append(line_list, i+1);
    else if (buffer->data[i] == '\0')
      break;
  }
  // Pick up the detritus at the end if necessary.
  if (buffer->data[buffer->size-2] != '\n')
    size_t_slist_append(line_list, buffer->size);

  buffer->num_lines = line_list->size;
  buffer->line_offsets = polymec_malloc(sizeof(size_t) * (buffer->num_lines+1));
  buffer->line_offsets[0] = 0;
  int i = 1;
  size_t line_no;
  size_t_slist_node_t* iter = NULL;
  while (size_t_slist_next(line_list, &iter, &line_no))
  {
    buffer->line_offsets[i] = line_no;
    ++i;
  }
  size_t_slist_free(line_list);
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
  buffer->size = (size_t)(end - start + 1);
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
  buffer->data[buffer->size-1] = '\0';

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

size_t text_buffer_size(text_buffer_t* buffer)
{
  return buffer->size - 1;
}

size_t text_buffer_num_lines(text_buffer_t* buffer)
{
  return buffer->num_lines;
}

bool text_buffer_next(text_buffer_t* buffer, int* pos, char** line, size_t* line_length)
{
  if (*pos >= (int)buffer->num_lines)
    return false;
  size_t offset = buffer->line_offsets[*pos];
  size_t next_offset = buffer->line_offsets[*pos+1];
  *line_length = next_offset - offset - 1;
  *line = &buffer->data[offset];
  ++(*pos);
  return true;
}

bool text_buffer_next_nonempty(text_buffer_t* buffer, int* pos, char** line, size_t* line_length)
{
  bool result = text_buffer_next(buffer, pos, line, line_length);
  while (result && (*line_length == 0))
    result = text_buffer_next(buffer, pos, line, line_length);
  return result;
}

char* text_buffer_to_string(text_buffer_t* buffer)
{
  return string_dup(buffer->data);
}

