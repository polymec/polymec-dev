// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <dlfcn.h>
#include "core/gfx.h"

typedef struct
{
  int page_index;  // Page index.
  int index; // Figure index within its page.
} gfx_figure_t;

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

static gfx_figure_t* gfx_current_figure(lua_State* L)
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

static gfx_figure_t* gfx_figure_new(int page_index, int fig_index)
{
  ASSERT(page_index > 0);
  ASSERT(fig_index > 0);

  gfx_figure_t* fig = polymec_malloc(sizeof(gfx_figure_t));
  fig->page_index = page_index;
  fig->index = fig_index;
  return fig;
}

gfx_figure_t* gfx_figure_new()
{
}

gfx_page_t* gfx_figure_page(gfx_figure_t* fig)
{
}

void gfx_figure_set_x_label(gfx_figure_t* fig, const char* label)
{
}

void gfx_figure_set_y_label(gfx_figure_t* fig, const char* label)
{
}

void gfx_figure_set_z_label(gfx_figure_t* fig, const char* label)
{
}

void gfx_figure_set_title(gfx_figure_t* fig, const char* title)
{
}

void gfx_figure_set_discrete_colors(gfx_figure_t* fig,
                                    int* rgba,
                                    size_t num_colors)
{
}

void gfx_figure_set_cont_colors(gfx_figure_t* fig,
                                int* rgba,
                                size_t num_points)
{
}

void gfx_figure_colorbar(gfx_figure_t* fig)
{
}

void gfx_figure_legend(gfx_figure_t* fig)
{
}

void gfx_figure_plot(gfx_figure_t* fig, 
                     real_t* x, 
                     real_t* y, 
                     size_t n,
                     char* label)
{
}

void gfx_figure_scatter(gfx_figure_t* fig,
                        real_t* x, 
                        real_t* y,
                        size_t n)
{
}

void gfx_figure_contour(gfx_figure_t* fig)
{
}

void gfx_figure_surface(gfx_figure_t* fig)
{
}

void gfx_figure_quiver(gfx_figure_t* fig)
{
}

void gfx_figure_image(gfx_figure_t* fig)
{
}

void gfx_figure_clear(gfx_figure_t* fig)
{
}

gfx_page_t* gfx_page_new(int num_rows, int num_cols)
{
}

int gfx_page_num_rows(gfx_page_t* page)
{
}

int gfx_page_num_columns(gfx_page_t* page)
{
}

bool gfx_page_next(gfx_page_t* page, int* pos, gfx_figure_t** figure)
{
}

gfx_figure_t* gfx_page_figure(gfx_page_t* page,
                              int row,
                              int column)
{
}

