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

#ifndef POLYMEC_HASH_FUNCTIONS_H
#define POLYMEC_HASH_FUNCTIONS_H

#include <string.h>

// A solid string hash function by Dan Bernstein.
static inline int djb2_hash(unsigned char* str, int len)
{
  int hash = 5381;

  for (int i = 0; i < len; ++i)
  {
    int c = str[i];
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */
  }

  return hash;
}

// A slightly improved version of djb2 using xor.
static inline int djb2_xor_hash(unsigned char* str, int len)
{
  int hash = 5381;

  for (int i = 0; i < len; ++i)
  {
    int c = str[i];
    hash = ((hash << 5) + hash) ^ c;
  }

  return hash;
}

static inline int int_hash(int i)
{
  return djb2_xor_hash((unsigned char*)&i, sizeof(int));
}

static inline int int_pair_hash(int* i)
{
  return djb2_xor_hash((unsigned char*)i, 2*sizeof(int));
}

static inline int string_hash(char* str)
{
  return djb2_xor_hash((unsigned char*)str, strlen(str));
}

static inline int pointer_hash(void* p)
{
  return djb2_xor_hash((unsigned char*)p, sizeof(void*));
}

#endif
