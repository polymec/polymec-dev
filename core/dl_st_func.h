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

#ifndef POLYMEC_DL_ST_FUNC_H
#define POLYMEC_DL_ST_FUNC_H

#include "st_func.h"

// Construct a space-time function from a dynamically-loaded module
// that is either compiled in-place or loaded into Polymec's runtime database.
st_func_t* dl_st_func_new(const char* name);

// Sets the C compiler to use to build dynamically-loaded functions.
void dl_st_func_set_compiler(const char* cc,
                             const char* cflags);

// Sets the directory in which shared objects are stored.
void dl_st_func_set_so_dir(const char* path);

// Registers the dynamically-loaded function with the given name with 
// Polymec, building its shared object using the given source code.
void dl_st_func_register(const char* name, const char* source_code);

#endif

