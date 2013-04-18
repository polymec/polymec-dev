// This implements polymesher's capability for writing meshes that can 
// be used by TOUGH2 and TOUGH+.

#include <string.h>
#include "core/interpreter.h"

#ifdef __cplusplus
extern "C" {
#endif

// Lua stuff.
#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

// write_tough_mesh(args) -- This function writes a given mesh to a file 
// on disk. Arguments (passed in a table according to Chapter 5.3 of the 
// Lua reference manual) are:
//
// filename -> name of the file to write (1 file only)
// format -> 'T2', 'T+' (mesh format to use)
// mesh -> mesh object 
// inactive_tag -> the tag within the mesh object that denotes inactive elements.
int write_tough_mesh(lua_State* lua)
{
  // Check the arguments.
  int num_args = lua_gettop(lua);
  if ((num_args != 1) || !lua_istable(lua, 1))
  {
    lua_pushstring(lua, "write_tough_mesh: invalid arguments. Usage:\n"
                        "write_tough_mesh{filename [= 'MESH'], format [= 'T2'/'T+'], mesh [= mesh], inactive_tag [= 'inactive']}).");
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Get the argument(s).
  mesh_t* mesh = NULL;
  lua_getfield(lua, 1, "mesh"); // required!
  if (!lua_isnoneornil(lua, 2))
    mesh = lua_tomesh(lua, 2);
  else
  {
    lua_pushstring(lua, "write_tough_mesh: mesh argument is required!");
    lua_error(lua);
    return LUA_ERRRUN;
  }
  lua_pop(lua, 1);

  char* filename = NULL;
  lua_getfield(lua, 1, "filename");
  if (!lua_isnoneornil(lua, 2))
    filename = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  char* format = NULL;
  lua_getfield(lua, 1, "format");
  if (!lua_isnoneornil(lua, 2))
    format = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  char* inactive_tag = NULL;
  lua_getfield(lua, 1, "inactive_tag");
  if (!lua_isnoneornil(lua, 2))
    inactive_tag = strdup(lua_tostring(lua, 2));
  lua_pop(lua, 1);

  // Provide defaults.
  if (filename == NULL)
    filename = strdup("MESH");
  if (format == NULL)
    format = strdup("T2");

  // Check our arguments.
  if ((strcasecmp(format, "t2") != 0) && (strcasecmp(format, "t+") != 0))
  {
    char err[128];
    snprintf(err, 128, "write_tough_mesh: unrecognized format: '%s'", format);
    lua_error(lua);
    return LUA_ERRRUN;
  }

  // Write the mesh to a file.
  // FIXME

  return 1;
}

#ifdef __cplusplus
}
#endif


