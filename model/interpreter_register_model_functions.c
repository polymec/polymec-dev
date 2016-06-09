// Copyright (c) 2012-2016, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "model/interpreter.h"
#include "model/periodic_bc.h"

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// Type code for periodic BCs.
extern int periodic_bc_type_code;

// Creates a periodic boundary condition from a pair of tags.
static int periodic_bc(lua_State* lua)
{
  // Check the argument.
  int num_args = lua_gettop(lua);
  if (num_args != 2)
    return luaL_error(lua, "Arguments must be 2 boundary mesh (face) tags.");

  for (int i = 1; i <= 3; ++i)
  {
    if (!lua_isstring(lua, i))
      return luaL_error(lua, "Argument %d must be a face tag.", i);
  }

  const char* tag1 = lua_tostring(lua, 1);
  const char* tag2 = lua_tostring(lua, 2);
  periodic_bc_t* bc = periodic_bc_new(tag1, tag2);
  lua_pushuserdefined(lua, bc, periodic_bc_type_code, NULL);
  return 1;
}

static docstring_t* periodic_bc_doc()
{
  return docstring_from_string("periodic_bc(face_tag1, face_tag2) - \n"
                               "  Creates a periodic boundary condition that aliases the sets of faces\n"
                               "  identified by the two face tags in a mesh.");
}

void interpreter_register_model_functions(interpreter_t* interp)
{
  periodic_bc_type_code = interpreter_new_user_defined_type_code(interp);
  interpreter_register_function(interp, "periodic_bc", periodic_bc, periodic_bc_doc());
}
