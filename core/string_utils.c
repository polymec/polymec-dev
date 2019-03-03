// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <string.h>
#include <ctype.h>
#include "core/string_utils.h"
#include "core/slist.h"

char* string_dup(const char* s)
{
  int len = (int)strlen(s);
  char* copy = polymec_malloc((len+1)*sizeof(char));
  strcpy(copy, s);
  return copy;
}

char* string_ndup(const char* s, size_t n)
{
  size_t len = MIN(n, strlen(s));
  char* copy = polymec_malloc((len+1)*sizeof(char));
  strncpy(copy, s, n);
  copy[n] = '\0';
  return copy;
}

int string_casecmp(const char* s1, const char* s2)
{
  int l1 = (int)strlen(s1), l2 = (int)strlen(s2);
  int lmin = MIN(l1, l2);
  for (int i = 0; i < lmin; ++i)
  {
    char c1 = (char)tolower(s1[i]), c2 = (char)tolower(s2[i]);
    if (c1 < c2)
      return -1;
    else if (c1 > c2)
      return 1;
  }
  if (l1 < l2)
    return -1;
  else if (l1 > l2)
    return 1;
  else
    return 0;
}

int string_ncasecmp(const char* s1, const char* s2, size_t n)
{
  int l1 = (int)strlen(s1), l2 = (int)strlen(s2);
  int lmin = MIN(l1, l2);
  if (lmin < n)
    return string_casecmp(s1, s2);
  for (int i = 0; i < n; ++i)
  {
    char c1 = (char)tolower(s1[i]), c2 = (char)tolower(s2[i]);
    if (c1 < c2)
      return -1;
    else if (c1 > c2)
      return 1;
  }
  return 0;
}

void string_free(char* s)
{
  polymec_free(s);
}

void string_copy_from_raw(const char* raw_array, size_t n, char* dest)
{
  strncpy(dest, raw_array, n);
  dest[n] = '\0';
}

bool string_next_token(const char* s, const char* delimiter, int* pos, char** token, size_t* length)
{
  int slen = (int)strlen(s);
  if (*pos >= slen)
    return false;
  char* delim = strstr((const char*)&(s[*pos]), delimiter);
  if ((delim == NULL) && (*pos == 0))
    return false;
  *token = (char*)&(s[*pos]);
  if (delim == NULL)
  {
    *length = (size_t)(slen - *pos);
    *pos = slen;
  }
  else
  {
    *length = (size_t)(delim - s) - *pos;
    *pos = (int)(delim - s) + (int)strlen(delimiter);
  }
  return true;
}

int string_num_tokens(const char* s, const char* delimiter)
{
  int num_tokens = 0, pos = 0;
  size_t length;
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

  int i = 0, pos = 0;
  size_t length;
  char* token;
  char** strs = polymec_malloc(sizeof(char*) * (*num_substrings));
  while (string_next_token(s, delimiter, &pos, &token, &length))
    strs[i++] = string_ndup(token, length);
  return strs;
}

int string_trim(char* s)
{
  int slen = (int)strlen(s);
  if (slen == 0) return 0;
  int l = 0, r = slen - 1;
  while ((r > 0) && isspace(s[r])) --r;
  while ((l < r) && isspace(s[l])) ++l;
  s[r+1] = '\0';
  return l;
}

bool string_is_number(const char* s)
{
  if ((s == NULL) || (*s == '\0') || isspace(*s))
    return false;
  char* p;
  strtod(s, &p);
  return (*p == '\0');
}

bool string_is_integer(const char* s)
{
  if ((s == NULL) || (*s == '\0') || isspace(*s))
    return false;
  char* p;
  strtol(s, &p, 10);
  return (*p == '\0');
}

bool string_as_boolean(const char* s)
{
  return ((s != NULL) &&
          ((strcmp(s, "1") == 0) ||
           (string_casecmp(s, "true") == 0) ||
           (string_casecmp(s, "yes") == 0) ||
           (string_casecmp(s, "on") == 0)));
}

int string_find_in_list(const char* s,
                        const char** string_list,
                        bool case_sensitive)
{
  ASSERT(string_list != NULL);
  if (s == NULL)
    return -1;

  int i = 0;
  int (*cmp)(const char* s1, const char* s2) = case_sensitive ? strcmp : string_casecmp;
  while (string_list[i] != NULL)
  {
    if (cmp(s, string_list[i]) == 0)
      return i;
    ++i;
  }
  return -1;
}

bool string_contains(const char* s, const char* subs)
{
  return (strstr(s, subs) != NULL);
}

// This stuff is used for string_substitute, below.

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

char* string_substitute(const char* string, string_substitution_t substitutions[])
{
  if (string == NULL)
    return NULL;
  else if (substitutions == NULL)
    return string_dup(string);

  // How many tokens have we got?
  int num_tokens = 0;
  while (strcmp(substitutions[num_tokens++].value, END_OF_SUBSTITUTION.value));

  // What's the length of each token / value?
  int token_lengths[num_tokens], value_lengths[num_tokens];
  for (int i = 0; i < num_tokens; ++i)
  {
    token_lengths[i] = (int)strlen(substitutions[i].token);
    value_lengths[i] = (int)strlen(substitutions[i].value);
  }

  // Now traverse the given string and jot down where each token occurs
  // (and which one it is).
  int_slist_t* token_occ = int_slist_new();
  int_slist_t* token_which = int_slist_new();
  for (int i = 0; i < num_tokens; ++i)
  {
    const char* token = substitutions[i].token;
    if ((token == NULL) || ((int)strlen(token) == 0)) continue;
    char* c = (char*)string;
    while (c != NULL)
    {
      c = strstr((const char*)c, token);
      if (c != NULL)
      {
        int_slist_append(token_occ, (int)(c - (char*)string));
        int_slist_append(token_which, i);
        c += 1;
      }
    }
  }
  size_t num_occ = token_occ->size;

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
  int old_len = (int)strlen(string);
  int new_len = old_len;
  for (int i = 0; i < num_occ; ++i)
  {
    int token_index = subst_data[i].which;
    new_len += value_lengths[token_index] - token_lengths[token_index];
  }
  new_len += 1; // '\0'

  // Now create the new string.
  char* new_string = polymec_malloc(sizeof(char)*new_len);
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

