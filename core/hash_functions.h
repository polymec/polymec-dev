#ifndef ARBI_HASH_FUNCTIONS_H
#define ARBI_HASH_FUNCTIONS_H

#ifdef __cplusplus
extern "C" {
#endif

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
  return djb2_xor_hash((unsigned char*)&i, sizeof(i));
}

static inline int string_hash(char* str)
{
  return djb2_xor_hash((unsigned char*)str, strlen(str));
}

#ifdef __cplusplus
}
#endif

#endif
