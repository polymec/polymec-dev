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

#ifndef POLYMEC_CONSTANT_ST_FUNC_H
#define POLYMEC_CONSTANT_ST_FUNC_H

#include "core/st_func.h"

// Construct a constant space-time function given a number of components 
// and their values.
st_func_t* constant_st_func_new(int num_comp, double comp[]);

// Free of charge, we toss in the sp_func version.
sp_func_t* constant_sp_func_new(int num_comp, double comp[]);

#endif

