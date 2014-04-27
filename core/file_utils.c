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

#include "core/polymec.h"
#include <unistd.h>
#include "file_utils.h"

#define SEPARATOR '/'

void parse_path(const char *path, char *dirname, char *filename)
{
  int len = strlen(path);
  char* last_sep = strrchr(path, SEPARATOR);
  // Now last_sep points to the last separator in path.

  if (last_sep == NULL) // No separator found!
  {
    // dirname is '.'
    dirname[0] = '.';
    dirname[1] = '\0';

    // filename is path.
    strcpy(filename, path);
  }
  else
  {
    int index = last_sep - path;
    strncpy(dirname, path, index);
    dirname[index] = '\0';

    strncpy(filename, path + index + 1, len - index - 1);
    filename[len - index - 1] = '\0';
  }
} 

void join_paths(const char *dirname, const char *filename, char *path)
{
  // If the directory includes a separator at the end, we don't add another one. 
  if (dirname[strlen(dirname)-1] == SEPARATOR)
    snprintf(path, FILENAME_MAX, "%s%s", dirname, filename);
  else
    snprintf(path, FILENAME_MAX, "%s%c%s", dirname, SEPARATOR, filename);
}

FILE* make_temp_file(char* filename)
{
  int fd = mkstemp(filename);
  if (fd == -1)
    polymec_error("make_temp_file: No file could be created.");
  return fdopen(fd, "w");
}

FILE* make_temp_file_with_suffix(char* filename, const char* suffix)
{
  char template[FILENAME_MAX];
  snprintf(template, FILENAME_MAX, "%s%s", filename, suffix);
  int fd = mkstemps(template, strlen(suffix));
  if (fd == -1)
    polymec_error("make_temp_file_with_suffix: No file could be created.");
  return fdopen(fd, "w");
}

