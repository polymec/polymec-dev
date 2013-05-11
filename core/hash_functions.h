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

#endif
