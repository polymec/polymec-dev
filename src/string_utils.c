#include <string.h>
#include "string_utils.h"

#define SEPARATOR '/'

#ifdef __cplusplus
extern "C" {
#endif

int parse_path(const char *path, char *dirname, char *filename)
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
  return ARBI_SUCCESS;
} 

int join_paths(const char *dirname, const char *filename, char *path)
{
  snprintf(path, 1024, "%s%c%s", dirname, SEPARATOR, filename);
  return ARBI_SUCCESS;
}

#ifdef __cplusplus
}
#endif

