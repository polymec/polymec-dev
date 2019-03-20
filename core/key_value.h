// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_KEY_VALUE_H
#define POLYMEC_KEY_VALUE_H

#include "core/polymec.h"

/// \addtogroup core core
///@{

/// \def DEFINE_KEY_VALUE(key_value_name, key_type, value_type)
/// A key_value is a pair associating a key with a value. It is used by
/// map-like containers. Data members for a key-value pair `kv`:
/// * `kv->key` - The key in the pair
/// * `kv->value` - The value in the pair
///
/// \param key_value_name The name of the key-value pair type.
/// \param key_type The data type for the key in the pair.
/// \param value_type The data type for the value in the pair.
#define DEFINE_KEY_VALUE(key_value_name, key_type, value_type) \
typedef struct \
{ \
  key_type key; \
  value_type value; \
} key_value_name##_t; \

///@}

// Some prototypical key-value pairs.
DEFINE_KEY_VALUE(string_int_key_value, char*, int)
DEFINE_KEY_VALUE(string_real_key_value, char*, real_t)
DEFINE_KEY_VALUE(string_ptr_key_value, char*, void*)

DEFINE_KEY_VALUE(int_int_key_value, int, int)
DEFINE_KEY_VALUE(int_real_key_value, int, real_t)
DEFINE_KEY_VALUE(int_ptr_key_value, int, void*)

#endif
