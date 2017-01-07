// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_DIFFERENCE_SP_FUNC_H
#define POLYMEC_DIFFERENCE_SP_FUNC_H

#include "core/sp_func.h"

// This signed distance function represents the difference of two 
// given surfaces represented by signed distance functions.
sp_func_t* difference_sp_func_new(sp_func_t* surface1, sp_func_t* surface2);

#endif

