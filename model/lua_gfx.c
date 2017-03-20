// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <dlfcn.h>
#include "core/lua_core.h"

#include "lua.h"
#include "lualib.h"
#include "lauxlib.h"

typedef struct
{
  int page_index;  // Page index.
  int index; // Figure index within its page.
} gfx_fig_t;

typedef struct
{
  int index;
  int nrows, ncols;
  ptr_array_t* figures;
  size_t current_fig;
} gfx_page_t;

// Global graphics context.
typedef struct
{
  lua_State* L;
  char path[FILENAME_MAX+1];

  // PlPlot functions.
  bool loaded;
  void* plplot;
  void (*plinit)();
  void (*plend)();
  void (*pladv)(int);
  void (*plssub)(int, int);
  void (*plstring)(int, double*, double*, const char*);

  // Pages.
  ptr_array_t* pages;
  size_t current_page;
} gfx_t;

static gfx_t* _gfx = NULL;

// Use this to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr, fail_label) \
  { \
    void* ptr = dlsym(dylib, symbol_name); \
    if (ptr == NULL) \
    { \
      log_urgent("%s: unable to find %s in dynamic library.", __func__, symbol_name); \
      goto fail_label; \
    } \
    *((void**)&(function_ptr)) = ptr; \
  } 

static void gfx_load()
{
  ASSERT(_gfx != NULL);

  // Try to discover the path to plplot.
  lua_getglobal(_gfx->L, "gfx");
  lua_getfield(_gfx->L, -1, "plplot_path");
  char plplot_path[FILENAME_MAX+1];
  if (!lua_isstring(_gfx->L, -1))
    strncpy(plplot_path, "INVALID_PLPLOT_PATH", FILENAME_MAX);
  else
    strncpy(plplot_path, lua_tostring(_gfx->L, -1), FILENAME_MAX);

  // Try to find the library.
  void* plplot = dlopen(plplot_path, RTLD_NOW);
#define FETCH_PLPLOT_SYMBOL(symbol_name) \
  FETCH_SYMBOL(plplot, #symbol_name, _gfx->symbol_name, failure);
  FETCH_PLPLOT_SYMBOL(plinit);
  FETCH_PLPLOT_SYMBOL(plend);
  FETCH_PLPLOT_SYMBOL(pladv);
  FETCH_PLPLOT_SYMBOL(plssub);
  FETCH_PLPLOT_SYMBOL(plstring);
#undef FETCH_PLPLOT_SYMBOL

  _gfx->plplot = plplot;
  _gfx->pages = ptr_array_new();
  strncpy(_gfx->path, plplot_path, FILENAME_MAX);
  _gfx->loaded = true;
  return;

failure:
  dlclose(plplot);
  log_info("Could not load plplot library from %s.", plplot_path);
  _gfx->loaded = false;
}

static void gfx_finalize()
{
  ptr_array_free(_gfx->pages);
  polymec_free(_gfx);
  _gfx = NULL;
}

static gfx_t* gfx_instance(lua_State* L)
{
  if (_gfx == NULL)
  {
    _gfx = polymec_malloc(sizeof(gfx_t));
    _gfx->L = L;
    _gfx->loaded = false;
    polymec_atexit(gfx_finalize);
  }

  // Keep trying to load.
  if (!_gfx->loaded)
    gfx_load();

  return _gfx;
}

static gfx_fig_t* gfx_current_figure(lua_State* L)
{
  gfx_t* gfx = gfx_instance(L);
  gfx_page_t* page = gfx->pages->data[gfx->current_page];
  return page->figures->data[page->current_fig];
}

static void gfx_page_free(gfx_page_t* page)
{
  ptr_array_free(page->figures);
  polymec_free(page);
}

static gfx_fig_t* gfx_figure_new(int page_index, int fig_index)
{
  ASSERT(page_index > 0);
  ASSERT(fig_index > 0);

  gfx_fig_t* fig = polymec_malloc(sizeof(gfx_fig_t));
  fig->page_index = page_index;
  fig->index = fig_index;
  return fig;
}

static gfx_page_t* gfx_page_new(lua_State* L, int nrows, int ncols)
{
  ASSERT(nrows > 0);
  ASSERT(ncols > 0);

  gfx_t* gfx = gfx_instance(L);
  gfx_page_t* page = polymec_malloc(sizeof(gfx_page_t));
  page->index = (int)gfx->pages->size+1;
  page->nrows = nrows;
  page->ncols = ncols;
  page->figures = ptr_array_new();
  for (int i = 0; i < nrows*ncols; ++i)
  {
    gfx_fig_t* fig = gfx_figure_new(page->index, i);
    ptr_array_append_with_dtor(page->figures, fig, polymec_free);
  }

  page->current_fig = 0;
  return page;
}

static int gfx_page_figure(lua_State* L)
{
  gfx_page_t* page = luaL_checkudata(L, 1, "gfx.page");

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

  if ((row < 1) || (row > page->nrows))
  {
    if (page->nrows > 1)
      return luaL_error(L, "row must be between 1 and %d.", page->nrows);
    else
      return luaL_error(L, "row must be between 1.");
  }
  if ((col < 1) || (col > page->ncols))
  {
    if (page->ncols > 1)
      return luaL_error(L, "columm must be between 1 and %d.", page->ncols);
    else
      return luaL_error(L, "column must be between 1.");
  }

  int fig_index = row*page->ncols + col;
  gfx_fig_t* fig = page->figures->data[fig_index];
  lua_push_object(L, "gfx.figure", fig, NULL);
  return 1;
}

static int gfx_page_tostring(lua_State* L)
{
  gfx_page_t* page = luaL_checkudata(L, 1, "gfx.page");
  lua_pushfstring(L, "gfx.page (%i figures)", (int)(page->figures->size));
  return 1;
}

static lua_class_method gfx_page_methods[] = {
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

static lua_class_method gfx_fig_methods[] = {
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

  gfx_t* gfx = gfx_instance(L);
  gfx_page_t* page = gfx_page_new(L, nrows, ncols);
  ptr_array_append_with_dtor(gfx->pages, page, DTOR(gfx_page_free));
  lua_push_object(L, "gfx.page", page, NULL);
  return 1;
}

static int gfx_plot(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_scatter(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_contour(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_surface(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_quiver(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_surface3d(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_image(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static int gfx_clear(lua_State* L)
{
  gfx_fig_t* fig = gfx_current_figure(L);
  if (fig != NULL)
  {
    lua_push_object(L, "gfx.figure", fig, NULL);
    return 1;
  }
  else
    return 0;
}

static lua_module_function gfx_functions[] = {
  {"page", gfx_page},
  {"plot", gfx_plot},
  {"scatter", gfx_scatter},
  {"contour", gfx_contour},
  {"surface", gfx_surface},
  {"quiver", gfx_quiver},
  {"surface3d", gfx_surface3d},
  {"image", gfx_image},
  {"clear", gfx_clear},
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

  // Add a candidate path to the gfx table.
  lua_pushstring(L, "/usr/local/lib/libplpath.so");
  lua_setfield(L, -2, "plplot_path");
}

