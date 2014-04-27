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

#include <dirent.h>
#include <sys/stat.h>
#include "core/polymec.h"
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

// Global temporary directory for the present process.
static char* polymec_temp_dir = NULL;

// This helper removes the present process's temporary directory on exit.
static void remove_polymec_temp_dir()
{
  ASSERT(polymec_temp_dir != NULL);

  log_debug("Deleting temporary directory '%s'.", polymec_temp_dir);
  remove_dir(polymec_temp_dir);
  free(polymec_temp_dir);
}

// This helper creates the temporary directory for the present polymec 
// process and then fills temp_dir with the name.
static void make_polymec_temp_dir(char* temp_dir)
{
  if (polymec_temp_dir == NULL)
  {
    char* tmpdir = getenv("TMPDIR");
    int pid = getpid();
    if (tmpdir != NULL)
      snprintf(temp_dir, FILENAME_MAX, "%s/polymec-%d", tmpdir, pid);
    else
      snprintf(temp_dir, FILENAME_MAX, "/tmp/polymec-%d", pid);
    polymec_temp_dir = string_dup(temp_dir);
    mkdir(temp_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    polymec_atexit(remove_polymec_temp_dir);
  }
  else
    strncpy(temp_dir, polymec_temp_dir, FILENAME_MAX);
}

FILE* make_temp_file(const char* file_template, char* filename)
{
  // Construct polymec's temporary directory/name.
  char temp_dir[FILENAME_MAX];
  make_polymec_temp_dir(temp_dir);

  // Construct the full temporary file path.
  char full_template[FILENAME_MAX];
  snprintf(full_template, FILENAME_MAX, "%s/%s", temp_dir, file_template);

  // Make the thing.
  int fd = mkstemp(full_template);
  if (fd == -1)
    return NULL;
  else
    strncpy(filename, full_template, FILENAME_MAX);
  return fdopen(fd, "w");
}

bool make_temp_dir(const char* dir_template, char* dirname)
{
  // Construct polymec's temporary directory/name.
  char temp_dir[FILENAME_MAX];
  make_polymec_temp_dir(temp_dir);

  // Construct the full temporary directory path.
  char full_template[FILENAME_MAX];
  snprintf(full_template, FILENAME_MAX, "%s/%s", temp_dir, dir_template);

  // Make the thing.
  char* dir = mkdtemp(full_template);
  if (dir == NULL)
    polymec_error("make_temp_dir: No directory could be created.");
  else
    strncpy(dirname, full_template, FILENAME_MAX);
  return (dir != NULL);
}

int remove_dir(const char* path)
{
  DIR* d = opendir(path);
  int path_len = strlen(path);
  int r = -1;

  if (d != NULL)
  {
    struct dirent *p;

    r = 0;
    while (!r && (p = readdir(d)))
    {
      // Skip . and .. entries.
      if (!strcmp(p->d_name, ".") || !strcmp(p->d_name, ".."))
        continue;

      int len = path_len + strlen(p->d_name) + 2; 
      char buf[FILENAME_MAX];
      struct stat statbuf;

      snprintf(buf, len, "%s/%s", path, p->d_name);

      int r2 = -1;
      if (!stat(buf, &statbuf))
      {
        if (S_ISDIR(statbuf.st_mode))
          r2 = remove_dir(buf);
        else
          r2 = remove(buf);
      }

      r = r2;
    }

    closedir(d);
  }

  // Don't forget to remove the directory itself.
  if (!r)
   r = rmdir(path);

  return r;
}

