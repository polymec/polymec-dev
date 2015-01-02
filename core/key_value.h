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
DEFINE_KEY_VALUE(string_real_key_value, char*, real_t)
DEFINE_KEY_VALUE(string_ptr_key_value, char*, void*)

DEFINE_KEY_VALUE(int_int_key_value, int, int)
DEFINE_KEY_VALUE(int_real_key_value, int, real_t)
DEFINE_KEY_VALUE(int_ptr_key_value, int, void*)

#endif
