// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/lua_core.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

typedef struct
{
  int hello;
} gfx_page_t;

typedef struct
{
  int page_index; 
} gfx_fig_t;

static int gfx_page_figure(lua_State* L)
{
  return 0;
}

static int gfx_page_tostring(lua_State* L)
{
  return 0;
}

static lua_class_method[] gfx_page_methods = {
  {"figure", gfx_page_figure},
  {"__tostring", gfx_page_tostring},
  {NULL, NULL}
};

static int gfx_fig_x_axis(lua_State* L)
{
  return 0;
}

static int gfx_fig_y_axis(lua_State* L)
{
  return 0;
}

static int gfx_fig_title(lua_State* L)
{
  return 0;
}

static int gfx_fig_colorbar(lua_State* L)
{
  return 0;
}

static int gfx_fig_legend(lua_State* L)
{
  return 0;
}

static int gfx_fig_tostring(lua_State* L)
{
  return 0;
}

static lua_class_method[] gfx_fig_methods = {
  {"x_axis", gfx_fig_x_axis},
  {"y_axis", gfx_fig_y_axis},
  {"title", gfx_fig_title},
  {"colorbar", gfx_fig_colorbar},
  {"legend", gfx_fig_legend},
  {"__tostring", gfx_fig_tostring},
  {NULL, NULL}
};

static int gfx_page(lua_State* L)
{
  int num_args = lua_gettop(L);
  if ((num_args != 1) && (num_args != 2))
    return luaL_error(L, "Arguments must be a 2-tuple containing (nrows, ncols).");
  if ((num_args == 1) && !lua_istable(L, 1))
    return luaL_error(L, "Argument 1 must be a 2-tuple containing (nrows, ncols).");
  if ((num_args == 2) && (!lua_isinteger(L, 1) || !lua_isinteger(L, 2))
    return luaL_error(L, "Arguments must be integers.");

  // FIXME: Parse nrows, ncols.
  int nrows = 1, ncols = 1;

  gfx_page_t* page = gfx_page_new(nrows, ncols);
  lua_push_gfx_page(L, page);
  return 1;
}

static int gfx_plot(lua_State* L)
{
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_scatter(lua_State* L)
{
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_contour(lua_State* L)
{
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_surface(lua_State* L)
{
  lua_push_gfx_figure(L, fig);
  return 1;
}

static int gfx_quiver(lua_State* L)
{
  lua_push_gfx_quiver(L, fig);
  return 1;
}

static int gfx_surface3d(lua_State* L)
{
  lua_push_gfx_surface3d(L, fig);
  return 1;
}

static int gfx_image(lua_State* L)
{
  lua_push_gfx_image(L, fig);
  return 1;
}

static lua_module_function[] gfx_functions = {
  {"page", gfx_page},
  {"plot", gfx_plot},
  {"scatter", gfx_scatter},
  {"contour", gfx_contour},
  {"surface", gfx_surface},
  {"quiver", gfx_quiver},
  {"surface3d", gfx_surface3d},
  {"image", gfx_image},
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

  // Add the global page to the table.
//  lua_pushstring(L, "byte");
//  lua_setfield(L, -2, "byte");
}

