// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#ifndef POLYMEC_STRING_UTILS_H
#define POLYMEC_STRING_UTILS_H

#include <stdbool.h>

// Since strdup() is not standard C, we provide a surrogate here.
char* string_dup(const char* s);

// This is a version of string_dup() that only copies n characters, also 
// tacking an '\0' to the end.
char* string_ndup(const char* s, int n);

// This function copies at most n-1 characters from a raw array to the string 
// dest. The raw array is not assumed to be NULL-terminated, but if it is, 
// only those characters preceding the NULL character will be copied. The 
// destination string dest is assumed to have enough storage to hold at least 
// n characters including a NULL termination character. After the copy is 
// performed, a NULL character will be written to the (n-1)th index of dest.
void string_copy_from_raw(const char* raw_array, int n, char* dest);

// This just calls polymec_free() to free a string, but can be used as a 
// convenience function for assigning destructors to strings in containers.
void string_free(char* s);

// This function allows one to traverse a string containing a number of 
// delimiters, reading off the tokens encountered in between the delimiters.
// It returns false if the string has been completely traversed, true otherwise.
// pos should be set to 0 to begin the traversal. The next token and its length
// are stored in *token and *length, respectively.
bool string_next_token(const char* s, const char* delimiter, int* pos, char** token, int* length);

// Returns the number of substrings (separated by the given delimiter)
// occur in the given string.
int string_num_tokens(const char* s, const char* delimiter);

// Split a string into substrings using the given delimiter, and return 
// an array of newly-allocated strings. num_substrings will contain the 
// length of this array. If the delimitor is not found, the whole string will 
// be a single substring.
char** string_split(const char* s, const char* delimiter, int* num_substrings);

// This function allows one to remove whitespace from the front and back of 
// the given NULL-terminated string s. It replaces the character after the 
// last non-whitespace character in s with '\0', and returns the index of the 
// first non-whitespace character in s.
int string_trim(char* s);

// Returns true if the given string is numeric, false if not.
bool string_is_number(const char* s);

// Returns true if the given string contains the given substring, 
// false if not.
bool string_contains(const char* s, const char* subs);

// Given a string and a NULL-terminated list of token-value pairs, 
// returns a newly-allocated string containing the original string with 
// all instances of tokens replaced with their substitution values. This 
// string must be freed with polymec_free() or string_free().
typedef struct { char* token; char* value; } string_subst_t;
static const string_subst_t END_OF_SUBST = {(char*)"END_TOKEN", (char*)"END_VALUE"};
char* string_subst(const char* string, string_subst_t substitutions[]);

#endif
