// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <string.h>
#include "core/string_utils.h"
#include "core/slist.h"

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
  snprintf(path, 1024, "%s%c%s", dirname, SEPARATOR, filename);
}

// This stuff is used for string_subst, below.

typedef struct
{
  int occ_index;
  int which;
} string_subst_data_t;

static int string_subst_data_cmp(const void* left, const void* right)
{
  string_subst_data_t* l = (string_subst_data_t*)left;
  string_subst_data_t* r = (string_subst_data_t*)right;
  return (l->occ_index < r->occ_index) ? -1 : (l->occ_index == r->occ_index) ? 0 : 1;
}

char* string_subst(const char* string, string_subst_t substitutions[])
{
  if (string == NULL) 
    return NULL;
  else if (substitutions == NULL)
    return strdup(string);

  // How many tokens have we got?
  int num_tokens = 0;
  while (strcmp(substitutions[num_tokens++].value, END_OF_SUBST.value));

  // What's the length of each token / value?
  int token_lengths[num_tokens], value_lengths[num_tokens];
  for (int i = 0; i < num_tokens; ++i)
  {
    token_lengths[i] = strlen(substitutions[i].token);
    value_lengths[i] = strlen(substitutions[i].value);
  }

  // Now traverse the given string and jot down where each token occurs 
  // (and which one it is).
  int_slist_t* token_occ = int_slist_new();
  int_slist_t* token_which = int_slist_new();
  for (int i = 0; i < num_tokens; ++i)
  {
    const char* token = substitutions[i].token;
    if ((token == NULL) || (strlen(token) == 0)) continue;
    char* c = (char*)string;
    while (c != NULL)
    {
      c = strstr((const char*)c, token);
      if (c != NULL)
      {
        int_slist_append(token_occ, c - (char*)string);
        int_slist_append(token_which, i);
        c += 1;
      }
    }
  }

  // Copy this information into an array and sort it.
  int num_occ = token_occ->size;
  string_subst_data_t subst_data[num_occ];
  {
    int_slist_node_t* occ = token_occ->front;
    int_slist_node_t* which = token_which->front;
    int i = 0;
    while (occ != NULL)
    {
      subst_data[i].occ_index = occ->value;
      subst_data[i].which = which->value;
      ++i;
      occ = occ->next;
      which = which->next;
    }
  }
  int_slist_free(token_occ);
  int_slist_free(token_which);
  qsort(subst_data, num_occ, sizeof(string_subst_data_t), string_subst_data_cmp);

  // Compute the length of the new string.
  int old_len = strlen(string);
  int new_len = old_len;
  for (int i = 0; i < num_occ; ++i)
  {
    int token_index = subst_data[i].which;
    new_len += value_lengths[token_index] - token_lengths[token_index];
  }
  new_len += 1; // '\0'

  // Now create the new string.
  char* new_string = malloc(sizeof(char)*new_len);
  int read_index = 0, write_index = 0;
  for (int i = 0; i < num_occ; ++i)
  {
    // Copy the stuff between tokens from the old string.
    int occ_index = subst_data[i].occ_index;

    // Calculate the number of bytes to read / write.
    int read_len = occ_index; 
    if (i > 0)
    {
      int prev_occ = subst_data[i-1].occ_index;
      int prev_token_index = subst_data[i-1].which;
      int prev_token_len = token_lengths[prev_token_index];
      read_len -= (prev_occ + prev_token_len);
    }
    memcpy(&new_string[write_index], &string[read_index], read_len*sizeof(char));
    read_index += read_len;
    write_index += read_len;

    // Write the value instead of the token.
    int token_index = subst_data[i].which;
    const char* value = substitutions[token_index].value;
    memcpy(&new_string[write_index], value, value_lengths[token_index]*sizeof(char));
    read_index += token_lengths[token_index];
    write_index += value_lengths[token_index];
  }
  // Pick up the rest of the string.
  if (read_index < old_len)
  {
    memcpy(&new_string[write_index], &string[read_index], (old_len-read_index)*sizeof(char));
    write_index += old_len - read_index;
    read_index = old_len;
  }

  // Add the newline.
  ASSERT(write_index == new_len-1);
  new_string[write_index] = '\0';

  return new_string;
}

