// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_driver.h"
#include "core/lua_core.h"
#include "geometry/lua_geometry.h"
#include "model/lua_model.h"

static int register_modules(lua_State* L)
{
  lua_register_core_modules(L);
  lua_register_geometry_modules(L);
  lua_register_model_modules(L);
  return 0;
}

int main(int argc, char* argv[]) 
{
  return lua_driver(argc, argv, register_modules);
}
