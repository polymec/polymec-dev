#ifndef POLYMEC_KEY_VALUE_H
#define POLYMEC_KEY_VALUE_H

#include "core/polymec.h"

// A key_value is a pair associating a key with a value. It is used by 
// map-like containers.
// One defines a key-value pair using
// DEFINE_KEY_VALUE(key_value_name, key_type, value_type)

#define DEFINE_KEY_VALUE(key_value_name, key_type, value_type) \
typedef struct \
{ \
  key_type key; \
  value_type value; \
} key_value_name##_t; \

// Some prototypical key-value pairs.
DEFINE_KEY_VALUE(string_int_key_value, char*, int)
DEFINE_KEY_VALUE(string_double_key_value, char*, double)
DEFINE_KEY_VALUE(string_ptr_key_value, char*, void*)

DEFINE_KEY_VALUE(int_int_key_value, int, int)
DEFINE_KEY_VALUE(int_double_key_value, int, double)
DEFINE_KEY_VALUE(int_ptr_key_value, int, void*)

#endif
