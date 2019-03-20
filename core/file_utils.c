// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <dirent.h>
#include <unistd.h> // for rmdir().
#include "core/polymec.h"
#include "core/file_utils.h"
#include "core/string_utils.h"
#include "core/slist.h"

#define SEPARATOR '/'

// FIXME: These are not standard C.
extern int mkstemp(char *template);
extern FILE* fdopen(int fildes, const char *mode);
extern char *mkdtemp(char *template);

void parse_path(const char *path, char *dirname, char *filename)
{
  int len = (int)strlen(path);
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
    int index = (int)(last_sep - path);
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

bool create_directory(const char* dirname, mode_t mode)
{
  // Create every directory in the given path.
  int pos = 0;
  char path[FILENAME_MAX];
  char* delim = strstr((const char*)&(dirname[pos]), "/");
  while (delim != NULL)
  {
    pos = (int)(delim - dirname) + 1; // Include '/'
    snprintf(path, pos, "%s", dirname);
    DIR* dir = opendir(path);
    if (dir == NULL)
      mkdir(path, mode);
    else
      closedir(dir);
    delim = strstr((const char*)&(dirname[pos]), "/");
  }
  snprintf(path, FILENAME_MAX, "%s", dirname);
  return (mkdir(path, mode) == 0);
}

bool file_exists(const char* filename)
{
  FILE* file = fopen(filename, "r");
  bool result = (file != NULL);
  if (result)
    fclose(file);
  return result;
}

bool directory_exists(const char* dirname)
{
  DIR* dir = opendir(dirname);
  bool result = (dir != NULL);
  if (dir != NULL)
    closedir(dir);
  return result;
}

string_slist_t* files_within_directory(const char* dirname)
{
  DIR* dir = opendir(dirname);
  if (dir == NULL)
    polymec_error("Directory '%s' does not exist.", dirname);
  int path_len = (int)strlen(dirname);
  struct dirent* p;
  string_slist_t* files = string_slist_new();
  while ((p = readdir(dir)))
  {
    int len = path_len + (int)strlen(p->d_name) + 2;
    char entry_name[FILENAME_MAX];
    snprintf(entry_name, len, "%s/%s", dirname, p->d_name);
    struct stat statbuf;
    if (stat(entry_name, &statbuf) == 0)
    {
      if (!S_ISDIR(statbuf.st_mode))
        string_slist_append_with_dtor(files, string_dup(p->d_name), string_free);
    }
  }
  closedir(dir);
  return files;
}

string_slist_t* directories_within_directory(const char* dirname)
{
  DIR* dir = opendir(dirname);
  if (dir == NULL)
    polymec_error("Directory '%s' does not exist.", dirname);
  int path_len = (int)strlen(dirname);
  struct dirent* p;
  string_slist_t* dirs = string_slist_new();
  while ((p = readdir(dir)))
  {
    int len = path_len + (int)strlen(p->d_name) + 2;
    char entry_name[FILENAME_MAX];
    snprintf(entry_name, len, "%s/%s", dirname, p->d_name);
    struct stat statbuf;
    if (stat(entry_name, &statbuf) == 0)
    {
      if (S_ISDIR(statbuf.st_mode) &&
          (strcmp(p->d_name, ".") != 0) &&
          (strcmp(p->d_name, "..") != 0))
        string_slist_append_with_dtor(dirs, string_dup(p->d_name), string_free);
    }
  }
  closedir(dir);
  return dirs;
}

// Global temporary directory for the present process.
static char* polymec_temp_dir = NULL;

// This helper removes the present process's temporary directory on exit.
static void remove_polymec_temp_dir()
{
  ASSERT(polymec_temp_dir != NULL);

  log_debug("Deleting temporary directory '%s'.", polymec_temp_dir);
  remove_directory(polymec_temp_dir);
  polymec_free(polymec_temp_dir);
}

// This helper creates the temporary directory for the present polymec
// process and then fills temp_dir with the name.
static void make_polymec_temp_dir(char* temp_dir)
{
  if (polymec_temp_dir == NULL)
  {
    char* tmpdir = getenv("TMPDIR");
    int pid = getpid();
    if ((tmpdir != NULL) && (strlen(tmpdir) > 0))
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
  {
    log_urgent("make_temp_file: Couldn't create temporary file at %s.", full_template);
    return NULL;
  }
  else
    strncpy(filename, full_template, FILENAME_MAX);
  return fdopen(fd, "w");
}

bool make_temp_directory(const char* dir_template, char* dirname)
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

bool remove_directory(const char* path)
{
  DIR* d = opendir(path);
  int path_len = (int)strlen(path);
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

      int len = path_len + (int)strlen(p->d_name) + 2;
      char buf[FILENAME_MAX];
      struct stat statbuf;

      snprintf(buf, len, "%s/%s", path, p->d_name);

      int r2 = -1;
      if (!stat(buf, &statbuf))
      {
        if (S_ISDIR(statbuf.st_mode))
          r2 = remove_directory(buf);
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

  return (r == 0);
}

