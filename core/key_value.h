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
