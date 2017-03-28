// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"
#include "core/lua_gfx.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

static int gfx_page_fig(lua_State* L)
{
  gfx_page_t* page = lua_to_gfx_page(L, 1);

  int num_args = lua_gettop(L);
  if ((num_args != 2) && (num_args != 3))
    return luaL_error(L, "Argument must be a 2-tuple containing (row, column).");
  if ((num_args == 2) && !lua_istable(L, 1))
    return luaL_error(L, "Argument 1 must be a 2-tuple containing (row, column).");
  if ((num_args == 3) && (!lua_isinteger(L, 1) || !lua_isinteger(L, 2)))
    return luaL_error(L, "Arguments must be integers.");

  // Parse row, col.
  int row = 0, col = 0;
  if (num_args == 2)
  {
    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Argument must be a 2-tuple containing (row, column).");
    row = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_pushinteger(L, 2);
    lua_gettable(L, -2);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Argument must be a 2-tuple containing (row, column).");
    col = (int)lua_tointeger(L, -1);
  }
  else
  {
    row = (int)lua_tointeger(L, 1);
    col = (int)lua_tointeger(L, 2);
  }

  int num_rows = gfx_page_num_rows(page);
  if ((row < 1) || (row > num_rows))
  {
    if (num_rows > 1)
      return luaL_error(L, "row must be between 1 and %d.", num_rows);
    else
      return luaL_error(L, "row must be between 1.");
  }
  int num_cols = gfx_page_num_columns(page);
  if ((col < 1) || (col > num_cols))
  {
    if (num_cols > 1)
      return luaL_error(L, "columm must be between 1 and %d.", num_cols);
    else
      return luaL_error(L, "column must be between 1.");
  }

  gfx_figure_t* fig = gfx_page_figure(page, row, col);
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_page_tostring(lua_State* L)
{
  gfx_page_t* page = luaL_checkudata(L, 1, "gfx.page");
  lua_pushfstring(L, "gfx.page (%i x %i figures)", gfx_page_num_rows(page),
                  gfx_page_num_columns(page));
  return 1;
}

static lua_class_method gfx_page_methods[] = {
  {"figure", gfx_page_fig},
  {"__tostring", gfx_page_tostring},
  {NULL, NULL}
};

static int gfx_fig_x_axis(lua_State* L)
{
//  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  return 0;
}

static int gfx_fig_y_axis(lua_State* L)
{
//  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  return 0;
}

static int gfx_fig_z_axis(lua_State* L)
{
//  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  return 0;
}

static int gfx_fig_title(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  int num_args = lua_gettop(L);
  if ((num_args != 2) || !lua_isstring(L, 2))
    return luaL_error(L, "Argument must be a title for the figure.");

  const char* title = lua_tostring(L, 2);
  gfx_figure_set_title(fig, title);
  return 0;
}

static int gfx_fig_plot(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_plot(fig, NULL, NULL, 0, "glyphy", 0, "label");
  return 0;
}

static int gfx_fig_contour(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_contour(fig);
  return 0;
}

static int gfx_fig_surface(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_surface(fig);
  return 0;
}

static int gfx_fig_quiver(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_quiver(fig);
  return 0;
}

static int gfx_fig_image(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_image(fig);
  return 0;
}

static int gfx_fig_clear(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_clear(fig);
  return 0;
}

static int gfx_fig_colorbar(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_colorbar(fig, 0.0, 0.0);
  return 0;
}

static int gfx_fig_legend(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  gfx_figure_legend(fig, 0.0, 0.0);
  return 0;
}

static int gfx_fig_tostring(lua_State* L)
{
  gfx_figure_t* fig = lua_to_gfx_figure(L, 1);
  lua_pushfstring(L, "gfx.figure('%s')", gfx_figure_title(fig));
  return 1;
}

static lua_class_method gfx_fig_methods[] = {
  {"x_axis", gfx_fig_x_axis},
  {"y_axis", gfx_fig_y_axis},
  {"z_axis", gfx_fig_z_axis},
  {"title", gfx_fig_title},
  {"plot", gfx_fig_plot},
  {"contour", gfx_fig_contour},
  {"surface", gfx_fig_surface},
  {"quiver", gfx_fig_quiver},
  {"image", gfx_fig_image},
  {"clear", gfx_fig_clear},
  {"colorbar", gfx_fig_colorbar},
  {"legend", gfx_fig_legend},
  {"__tostring", gfx_fig_tostring},
  {NULL, NULL}
};

static int gfx_figure(lua_State* L)
{
  gfx_figure_t* fig = gfx_figure_new();
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_page(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && (num_args != 2))
    return luaL_error(L, "Arguments must be a 2-tuple containing (nrows, ncols).");
  if ((num_args == 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument must be a 2-tuple containing (nrows, ncols).");
  if ((num_args == 2) && (!lua_isinteger(L, 1) || !lua_isinteger(L, 2)))
    return luaL_error(L, "Arguments must be integers.");

  // Parse nrows, ncols.
  int nrows = 0, ncols = 0;
  if (num_args == 1)
  {
    lua_pushinteger(L, 1);
    lua_gettable(L, -2);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Argument must be a 2-tuple containing (nrows, ncols).");
    nrows = (int)lua_tointeger(L, -1);
    lua_pop(L, 1);
    lua_pushinteger(L, 2);
    lua_gettable(L, -2);
    if (!lua_isinteger(L, -1))
      return luaL_error(L, "Argument must be a 2-tuple containing (row, column).");
    ncols = (int)lua_tointeger(L, -1);
  }
  else
  {
    nrows = (int)lua_tointeger(L, 1);
    ncols = (int)lua_tointeger(L, 2);
  }

  gfx_page_t* page = gfx_page_new(nrows, ncols);
  lua_push_gfx_page(L, page);
  return 1;
}

static lua_module_function gfx_functions[] = {
  {"page", gfx_page},
  {"figure", gfx_figure},
  {NULL, NULL}
};

void lua_register_gfx(lua_State* L);
void lua_register_gfx(lua_State* L)
{
  // We define a gfx.page class representing a page containing multiple plots.
  lua_register_class(L, "gfx.page", NULL, gfx_page_methods);

  // We define a gfx.figure class representing a plot.
  lua_register_class(L, "gfx.figure", NULL, gfx_fig_methods);

  // Define the gfx module itself.
  lua_register_module(L, "gfx", gfx_functions);
}

void lua_push_gfx_figure(lua_State* L, gfx_figure_t* fig)
{
  lua_push_object(L, "gfx.figure", fig, NULL);
}

bool lua_is_gfx_figure(lua_State* L, int index)
{
  return lua_is_object(L, index, "gfx.figure");
}

gfx_figure_t* lua_to_gfx_figure(lua_State* L, int index)
{
  return (gfx_figure_t*)lua_to_object(L, index, "gfx.figure");
}

void lua_push_gfx_page(lua_State* L, gfx_page_t* page)
{
  lua_push_object(L, "gfx.page", page, NULL);
}

bool lua_is_gfx_page(lua_State* L, int index)
{
  return lua_is_object(L, index, "gfx.page");
}

gfx_page_t* lua_to_gfx_page(lua_State* L, int index)
{
  return (gfx_page_t*)lua_to_object(L, index, "gfx.page");
}

