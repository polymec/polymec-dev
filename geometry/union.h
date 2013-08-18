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

#ifndef POLYMEC_UNION_H
#define POLYMEC_UNION_H

#include "core/sp_func.h"

// This signed distance function represents the union of the set of 
// given surfaces represented by signed distance functions.
sp_func_t* union_new(sp_func_t** surfaces, int num_surfaces);

// This is a shorthand function that creates the union of two 
// surfaces.
sp_func_t* union_new2(sp_func_t* surface1, sp_func_t* surface2);

#endif

