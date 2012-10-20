#ifndef ARBI_HASH_FUNCTIONS_H
#define ARBI_HASH_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

// A solid string hash function by Dan Bernstein.
static inline int djb2_hash(char* str)
{
  int hash = 5381;
  int c;

  while ((c = *str++))
    hash = ((hash << 5) + hash) + c; /* hash * 33 + c */

  return hash;
}

// A slightly improved version of djb2 using xor.
static inline int djb2_xor_hash(char* str)
{
  int hash = 5381;
  int c;

  while ((c = *str++))
    hash = ((hash << 5) + hash) ^ c;

  return hash;
}

static inline int string_hash(char* str)
{
  return djb2_xor_hash(str);
}

#ifdef __cplusplus
}
#endif

#endif
