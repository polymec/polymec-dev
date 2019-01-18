-- This script generates the C code for the thirdparty module that exposes 
-- the metadata within 3rdparty/thirdparty.lua
--
-- Usage: lua generate_thirdparty_module.lua /path/to/thirdparty.lua /path/to/lua_thirdparty.c

require('io')

-- Read 3rdparty/thirdparty.lua
dofile(arg[1])

-- Open our generated C file for writing.
f = io.open(arg[2], 'w')

-- Generate boilerplate C code for the function declaration/definition.
f:write("// This file is generated by polymec's build system. Please don't touch it.\n\n")
f:write('#include "lua.h"\n#include "lualib.h"\n#include "lauxlib.h"\n\n')
f:write('int register_thirdparty_module(lua_State* L);\n')
f:write('int register_thirdparty_module(lua_State* L)\n{\n')
f:write('  // thirdparty module table\n')
f:write('  lua_newtable(L);\n\n')

-- This iterator function traverses a table in sorted key order.
function sorted_pairs(t)
  local a = {}
  for n in pairs(t) do table.insert(a, n) end
  table.sort(a)
  local i = 0
  local iter = function()
    i = i+1
    if a[i] == nil then return nil
    else return a[i], t[a[i]]
    end
  end
  return iter
end

-- For each of the modules, generate code.
libs = {}
for lib_name, lib_metadata in sorted_pairs(thirdparty) do
  f:write(string.format('  // %s\n', lib_name))
  f:write('  lua_newtable(L);\n')
  for key, value in sorted_pairs(lib_metadata) do
    if type(value) == 'string' then
      value = value:gsub('\n', '\\n'):gsub('"', '\\"')
      f:write(string.format('  lua_pushstring(L, "%s");\n', value))
    elseif type(value) == 'number' then
      f:write(string.format('  lua_pushinteger(L, %d);\n', value))
    elseif type(value) == 'table' then
      f:write('  lua_newtable(L);\n')
      for k, v in sorted_pairs(value) do
        if type(v) == 'string' then
          v= v:gsub('\n', '\\n\n\\'):gsub('"', '\\"')
          f:write(string.format('  lua_pushstring(L, "%s");\n', v))
        elseif type(v) == 'number' then
          f:write(string.format('  lua_pushinteger(L, %d);\n', v))
        end
        f:write(string.format('  lua_setfield(L, -2, "%s");\n', k))
      end
    end
    f:write(string.format('  lua_setfield(L, -2, "%s");\n', key))
  end
  f:write(string.format('  lua_setfield(L, -2, "%s");\n\n', lib_name))
end

-- Finish up and write the thing out.
f:write('  lua_setglobal(L, "thirdparty");\n')
f:write('  return 0;\n}\n')
f:close()
