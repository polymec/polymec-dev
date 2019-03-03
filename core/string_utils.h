// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_STRING_UTILS_H
#define POLYMEC_STRING_UTILS_H

#include <stdbool.h>

/// \addtogroup core core
///@{

/// Since strdup() is not standard C, we provide a surrogate here.
char* string_dup(const char* s);

/// This is a version of string_dup() that only copies n characters, also
/// tacking an '\0' to the end.
char* string_ndup(const char* s, size_t n);

///@{
/// Since strcasecmp and strncasecmp are not standard C, we provide surrogates
/// here.
int string_casecmp(const char* s1, const char* s2);
int string_ncasecmp(const char* s1, const char* s2, size_t n);
///@}

/// This function copies at most n characters from a raw array to the string
/// dest, which should be large enough to store n+1 characters. The raw array
/// is not assumed to be NULL-terminated, but if it is, only those characters
/// preceding the NULL character will be copied. After the copy is performed,
/// a NULL character will be written to dest[n].
void string_copy_from_raw(const char* raw_array, size_t n, char* dest);

/// This just calls polymec_free() to free a string, but can be used as a
/// convenience function for assigning destructors to strings in containers.
void string_free(char* s);

/// This function allows one to traverse a string containing a number of
/// delimiters, reading off the tokens encountered in between the delimiters.
/// \param [in] s The string containing delimiters and tokens.
/// \param [in] delimiter The delimiter that separates tokens in the string.
/// \param [in,out] pos Set pos to 0 to begin the traversal.
/// \param [out] token Stores a internal pointer to the next token in the traversal.
///                    This pointer is owned by the original string, so the caller must copy the
///                    data out of here to obtain a separate string.
/// \param [out] length Stores the length of the next token in the traversal.
/// \returns false if the string does not contain the given token, or if it has
/// been completely traversed. Otherwise returns true.
bool string_next_token(const char* s, const char* delimiter, int* pos, char** token, size_t* length);

/// Returns the number of substrings (separated by the given delimiter)
/// occur in the given string.
int string_num_tokens(const char* s, const char* delimiter);

/// Splits a string into substrings using the given delimiter, returning an array of substrings
/// found. If the delimitor is not found, the whole string will be a single substring and is
/// returned in an array of length 1.
/// \param [in] s The string to split.
/// \param [in] delimiter The delimiter separating the desired substrings.
/// \param [out] num_substrings stores the length of the returned array of strings.
/// \returns an array of newly-allocated strings created from the split.
char** string_split(const char* s, const char* delimiter, int* num_substrings);

/// This function allows one to remove whitespace from the front and back of
/// the given NULL-terminated string s. It replaces the character after the
/// last non-whitespace character in s with '\0', and returns the index of the
/// first non-whitespace character in s.
int string_trim(char* s);

/// Returns true if the given string is numeric, false if not.
bool string_is_number(const char* s);

/// Returns true if the given string can be interpreted as an integer,
/// false if not.
bool string_is_integer(const char* s);

/// Returns true if the given string is "1", "TRUE", "YES", "ON", or any
/// case-insensitive version thereof, and false otherwise.
bool string_as_boolean(const char* s);

/// Returns the index of the entry in the (NULL-terminated) string list that
/// matches the given string s, or -1 if s does not appear in string_list.
/// The case_sensitive flag determines whether the case has to match.
int string_find_in_list(const char* s,
                        const char** string_list,
                        bool case_sensitive);

/// Returns true if the given string contains the given substring,
/// false if not.
bool string_contains(const char* s, const char* subs);

/// \struct string_substitution
/// A token/value pair for performing string substitutions using
/// \ref string_substitute.
typedef struct { char* token; char* value; } string_substitution_t;

/// This indicates the terminus of a list of string substitutions to use
/// with \ref string_substitute.
static const string_substitution_t END_OF_SUBSTITUTION = {(char*)"END_TOKEN", (char*)"END_VALUE"};

/// Given a string and a NULL-terminated list of token-value pairs,
/// returns a newly-allocated string containing the original string with
/// all instances of tokens replaced with their substitution values. This
/// string must be freed with polymec_free() or string_free().
char* string_substitute(const char* string, string_substitution_t substitutions[]);

///@}

#endif
