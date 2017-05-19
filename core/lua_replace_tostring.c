// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "array.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static void append_number_value(lua_State* L, 
                                int index,
                                char_array_t* str)
{
  char val_str[33];
  if (lua_isinteger(L, index))
  {
    int val = (int)(lua_tointeger(L, index));
    snprintf(val_str, 32, "%d", val);
  }
  else
  {
    double val = lua_tonumber(L, index);
    snprintf(val_str, 32, "%g", val);
  }
  size_t val_len = strlen(val_str);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val_str, sizeof(char)*val_len);
}

static void append_string_value(lua_State* L, 
                                int index,
                                char_array_t* str)
{
  const char* val = lua_tostring(L, index);
  size_t val_len = strlen(val);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val, sizeof(char)*val_len);
}

static bool has_tostring(lua_State* L,
                         int index)
{
  // Does the data type have a metatable with a __tostring metamethod?
  bool has_tostring_mm = false;
  bool has_metatable = (lua_getmetatable(L, index) != 0);
  if (has_metatable)
  {
    lua_getfield(L, -1, "__tostring");
    if (lua_isfunction(L, -1))
    {
      has_metatable = true;
      lua_pop(L, 1);
    }
    lua_pop(L, 1); // Pop the metatable.
  }
  return has_tostring_mm;
}

static void append_value_using_tostring(lua_State* L, 
                                        int index, 
                                        char_array_t* str)
{
  lua_getmetatable(L, index);
  lua_getfield(L, -1, "__tostring");
  lua_pushvalue(L, index);
  lua_call(L, 1, 1);
  const char* val_str = lua_tostring(L, -1);
  size_t val_len = strlen(val_str);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val_str, sizeof(char)*val_len);
  lua_pop(L, 2);
}

static void append_function_value(lua_State* L, 
                                  int index, 
                                  char_array_t* str)
{
  const char* val_str;
  if (lua_iscfunction(L, index))
    val_str = "cfunction";
  else
    val_str = "function";
  size_t val_len = strlen(val_str);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val_str, sizeof(char)*val_len);
}

static void append_lightuserdata_value(lua_State* L, 
                                       int index, 
                                       char_array_t* str)
{
  void* ptr = lua_touserdata(L, index);
  char val_str[129];
  snprintf(val_str, 128, "lightuserdata (%p)", ptr);
  size_t val_len = strlen(val_str);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val_str, sizeof(char)*val_len);
}

static void append_userdata_value(lua_State* L, 
                                  int index, 
                                  char_array_t* str)
{
  void* ptr = lua_touserdata(L, index);
  char val_str[129];
  snprintf(val_str, 128, "userdata (%p)", ptr);
  size_t val_len = strlen(val_str);
  char_array_resize(str, str->size + val_len);
  memcpy(&str->data[str->size-val_len], val_str, sizeof(char)*val_len);
}

static void append_table_value(lua_State* L, 
                               int table_index, 
                               int indents,
                               char_array_t* str)
{
  // Get the f'reals table index.
  if (table_index < 0)
    table_index += lua_gettop(L) + 1;

  // Opening brace.
  for (int i = 0; i < indents; ++i)
    char_array_append(str, ' ');
  char_array_append(str, '{');

  // Pick through the table.
  int num_items = 0;
  bool show_as_list = true;
  lua_pushnil(L);
  while (lua_next(L, table_index))
  {
    ++num_items;

    // Key is at index -2, value is at -1.
    int ind = 0;
    if (lua_isinteger(L, -2))
      ind = (int)(lua_tointeger(L, -2));

    // We show this table as a list if all its keys are consecutive integers.
    if ((ind != num_items) && (show_as_list == true))
    {
      show_as_list = false;
      if (ind == 0)
        char_array_append(str, ',');
    }

    // Write out the key if we're not displaying in list mode.
    if (!show_as_list)
    {
      // Start a new line and indent.
      char_array_append(str, '\n');
      for (int i = 0; i < indents+1; ++i)
        char_array_append(str, ' ');

      // Write the value.
      if (lua_isnumber(L, -1))
        append_number_value(L, -1, str);
      else if (lua_isstring(L, -1))
        append_string_value(L, -1, str);
      else if (lua_isfunction(L, -1))
        append_function_value(L, -1, str);
      else if (lua_islightuserdata(L, -1))
        append_lightuserdata_value(L, -1, str);
      else if (has_tostring(L, -1))
        append_value_using_tostring(L, -1, str);
      else if (lua_isuserdata(L, -1)) // Userdata without metatable
        append_userdata_value(L, -1, str);
      else // Table without metatable
        append_table_value(L, -1, 0, str);

      char_array_append(str, ' ');
      char_array_append(str, '=');
      char_array_append(str, ' ');
    }
    else // Otherwise we have to make allowances for a preceding comma.
    {
      if (ind > 1)
      {
        char_array_append(str, ',');
        char_array_append(str, ' ');
      }
    }

    // Write the value.
    if (lua_isnumber(L, -1))
      append_number_value(L, -1, str);
    else if (lua_isstring(L, -1))
      append_string_value(L, -1, str);
    else if (lua_isfunction(L, -1))
      append_function_value(L, -1, str);
    else if (lua_islightuserdata(L, -1))
      append_lightuserdata_value(L, -1, str);
    else if (has_tostring(L, -1))
      append_value_using_tostring(L, -1, str);
    else if (lua_isuserdata(L, -1)) // Userdata without metatable
      append_userdata_value(L, -1, str);
    else // Table without metatable
      append_table_value(L, -1, 0, str);


    // Pop the value to fetch the next key.
    lua_pop(L, 1);
  }

  // Closing brace.
  if (!show_as_list)
  {
    char_array_append(str, '\n');
    for (int i = 0; i < indents; ++i)
      char_array_append(str, ' ');
  }
  char_array_append(str, '}');
}

// This is a fancy replacement for Lua's native tostring() function.
static int lua_tostring_replacement(lua_State* L)
{
  if (lua_istable(L, 1))
  {
    if (lua_getmetatable(L, 1) == 0)
    {
      // Here's a table with no metatable. Let's asplode its contents.
      char_array_t* str = char_array_new();
      append_table_value(L, -1, 0, str);
      char_array_append(str, '\0');

      // Push the string to the stack and get on with it.
      lua_pushstring(L, str->data); 
      char_array_free(str);
      return 1;
    }
    else
    {
      // Pop the metatable and proceed.
      lua_pop(L, 1);
    }
  }

  // Use the usual tostring() function. 
  lua_getglobal(L, "_tostring");
  lua_pushvalue(L, -2);
  lua_call(L, 1, 1);
  return 1;
}

void lua_replace_tostring(lua_State* L);
void lua_replace_tostring(lua_State* L)
{
  // Replace the tostring function (displays table contents).
  lua_getglobal(L, "tostring");
  lua_setglobal(L, "_tostring");
  lua_pushcfunction(L, lua_tostring_replacement);
  lua_setglobal(L, "tostring");
}
