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

#include <string.h>
#include <ctype.h>
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
  // If the directory includes a separator at the end, we don't add another one. 
  if (dirname[strlen(dirname)-1] == SEPARATOR)
    snprintf(path, 1024, "%s%s", dirname, filename);
  else
    snprintf(path, 1024, "%s%c%s", dirname, SEPARATOR, filename);
}

char* string_dup(const char* s)
{
  int len = strlen(s);
  char* copy = malloc((len+1)*sizeof(char));
  strcpy(copy, s);
  return copy;
}

char* string_ndup(const char* s, int n)
{
  int len = MIN(n, strlen(s));
  char* copy = malloc((len+1)*sizeof(char));
  strncpy(copy, s, n);
  copy[n] = '\0';
  return copy;
}

bool string_next_token(const char* s, const char* delimiter, int* pos, char** token, int* length)
{
  if (*pos >= strlen(s))
    return false;
  char* delim = strstr((const char*)&(s[*pos]), delimiter);
  *token = (char*)&(s[*pos]);
  if (delim == NULL)
  {
    *length = strlen(s) - *pos;
    *pos = strlen(s);
  }
  else
  {
    *length = (delim - s) - *pos;
    *pos = (delim - s) + strlen(delimiter);
  }
  return true;
}

int string_num_tokens(const char* s, const char* delimiter)
{
  int num_tokens = 0, pos = 0, length;
  char* token;
  while (string_next_token(s, delimiter, &pos, &token, &length))
    ++num_tokens;
  return num_tokens;
}

char** string_split(const char* s, const char* delimiter, int* num_substrings)
{
  *num_substrings = string_num_tokens(s, delimiter);
  if (*num_substrings == 0)
    return NULL;

  int i = 0, pos = 0, length;
  char* token;
  char** strs = malloc(sizeof(char*) * (*num_substrings));
  while (string_next_token(s, delimiter, &pos, &token, &length))
    strs[i++] = string_ndup(token, length);
  return strs;
}

bool string_is_number(const char* s)
{
  if ((s == NULL) || (*s == '\0') || isspace(*s))
    return 0;
  char* p;
  strtod(s, &p);
  return (*p == '\0');
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
    return string_dup(string);

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
  int num_occ = token_occ->size;

  // If there are no occurrences of the tokens, we simply copy 
  // the original string.
  if (num_occ == 0)
  {
    int_slist_free(token_occ);
    int_slist_free(token_which);
    return string_dup(string);
  }

  // Copy this information into an array and sort it.
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
//    read_index = old_len;
  }

  // Add the newline.
  ASSERT(write_index == new_len-1);
  new_string[write_index] = '\0';

  return new_string;
}

