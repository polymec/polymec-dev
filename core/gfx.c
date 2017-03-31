// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include <dlfcn.h>
#include "mpi.h"
#include "core/gfx.h"
#include "core/array.h"
#include "core/unordered_map.h"

// Global graphics context.
typedef struct
{
  bool loaded;

  // PlPlot library and functions.
  void* plplot;
  void (*plsdev)(const char*);
  void (*plinit)();
  void (*plend)();
  void (*plsstrm)(int);
  void (*plstar)(int, int);
  void (*pladv)(int);
  void (*pleop)();
  void (*plssub)(int, int);
  void (*plline)(int, double*, double*);
  void (*plstring)(int, double*, double*, const char*);
  void (*plsfont)(int, int, int);
  void (*plenv0)(double, double, double, double, int, int);
  void (*pllab)(const char*, const char*, const char*);
  void (*plcol0)(int);
  void (*plscmap0a)(int*, int*, int*, double*, int);
  void (*pllegend)(double*, double*, int, int, double, double, double, 
                   int, int, int, int, int, int, 
                   int*, double, double, double, double, 
                   int*, char**, int*, int*, double*, double*, 
                   int*, int*, double*, int*, double*, int*, char**);

  ptr_array_t* pages;
  string_ptr_unordered_map_t* colormaps;
  string_ptr_unordered_map_t* palettes;
} gfx_t;

static gfx_t* _gfx = NULL;

typedef struct
{
  char* label;
  char* glyph;
  int color;
} legend_item_t;

static legend_item_t* legend_item_new(const char* label, const char* glyph, int color)
{
  legend_item_t* item = polymec_malloc(sizeof(legend_item_t));
  item->label = string_dup(label);
  item->glyph = string_dup(glyph);
  item->color = color;
  return item;
}

static void legend_item_free(legend_item_t* item)
{
  string_free(item->label);
  string_free(item->glyph);
  polymec_free(item);
}

struct gfx_figure_t
{
  int index;
  gfx_page_t* page;

  char* title;
  char* x_label;
  char* y_label;
  char* z_label;
  real_t x_min, x_max, y_min, y_max;
  gfx_scale_t x_scale_type, y_scale_type;
  ptr_array_t* legend;
};

struct gfx_page_t
{
  int index;
  int num_rows, num_cols;
  ptr_array_t* figures;
};

// Use this to retrieve symbols from dynamically loaded libraries.
#define FETCH_SYMBOL(dylib, symbol_name, function_ptr, fail_label) \
  { \
    void* ptr = polymec_dlsym(dylib, "c_" symbol_name); \
    if (ptr == NULL) \
    { \
      log_urgent("%s: unable to find %s in dynamic library.", __func__, symbol_name); \
      goto fail_label; \
    } \
    *((void**)&(function_ptr)) = ptr; \
  }

static int _mpi_rank = -1;

static void gfx_load()
{
  ASSERT(_gfx != NULL);
  ASSERT(_gfx->loaded == false);

  if (_mpi_rank == 0)
  {
    // Try to load the library.
    void* plplot = polymec_dlopen("plplot");
    if (plplot == NULL)
      goto failure;

#define FETCH_PLPLOT_SYMBOL(symbol_name) \
  FETCH_SYMBOL(plplot, #symbol_name, _gfx->symbol_name, failure)
    FETCH_PLPLOT_SYMBOL(plsdev);
    FETCH_PLPLOT_SYMBOL(plinit);
    FETCH_PLPLOT_SYMBOL(plend);
    FETCH_PLPLOT_SYMBOL(plsstrm);
    FETCH_PLPLOT_SYMBOL(plstar);
    FETCH_PLPLOT_SYMBOL(pladv);
    FETCH_PLPLOT_SYMBOL(pleop);
    FETCH_PLPLOT_SYMBOL(plssub);
    FETCH_PLPLOT_SYMBOL(plline);
    FETCH_PLPLOT_SYMBOL(plstring);
    FETCH_PLPLOT_SYMBOL(plsfont);
    FETCH_PLPLOT_SYMBOL(plenv0);
    FETCH_PLPLOT_SYMBOL(pllab);
    FETCH_PLPLOT_SYMBOL(plcol0);
    FETCH_PLPLOT_SYMBOL(plscmap0a);
    FETCH_PLPLOT_SYMBOL(pllegend);
#undef FETCH_PLPLOT_SYMBOL

    _gfx->plplot = plplot;
#ifdef LINUX
    _gfx->plsdev("wxwidgets");
#else
    _gfx->plsdev("cairo");
#endif
    _gfx->plinit();
    _gfx->colormaps = string_ptr_unordered_map_new();
    _gfx->palettes = string_ptr_unordered_map_new();

    // Tell all our friends we've loaded.
    _gfx->loaded = true;
    MPI_Bcast(&_gfx->loaded, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  else
  {
    // Ask process 0 whether we've loaded.
    int loaded;
    MPI_Bcast(&loaded, 1, MPI_INT, 0, MPI_COMM_WORLD);
    _gfx->loaded = loaded;
    if (!_gfx->loaded)
      goto failure;
  }

  // Set a default font.
  gfx_font_t font = {.family = GFX_FONT_SERIF,
                     .style = GFX_FONT_UPRIGHT,
                     .weight = GFX_FONT_MEDIUM};
  gfx_set_font(font);
  return;

failure:
  {
    int zero = 0;
    MPI_Bcast(&zero, 1, MPI_INT, 0, MPI_COMM_WORLD);
  }
  if (_gfx->plplot != NULL)
    polymec_dlclose(_gfx->plplot);
  log_detail("gfx_load: Could not load plplot library.");
  _gfx->loaded = false;
}

static void gfx_finalize()
{
  ptr_array_free(_gfx->pages);
  if (_gfx->plplot != NULL)
  {
    string_ptr_unordered_map_free(_gfx->palettes);
    string_ptr_unordered_map_free(_gfx->colormaps);
    _gfx->plend();
    polymec_dlclose(_gfx->plplot);
  }
  polymec_free(_gfx);
  _gfx = NULL;
}

static gfx_t* gfx_instance()
{
  if (_gfx == NULL)
  {
    MPI_Comm_rank(MPI_COMM_WORLD, &_mpi_rank);
    _gfx = polymec_malloc(sizeof(gfx_t));
    _gfx->loaded = false;
    _gfx->plplot = NULL;
    _gfx->pages = ptr_array_new();
    polymec_atexit(gfx_finalize);
  }

  // Keep trying to load till we get it right!
  if (_gfx != NULL)
  {
    if (!_gfx->loaded)
      gfx_load();
  }

  return _gfx;
}

bool gfx_enabled()
{
  return ((gfx_instance() != NULL) && _gfx->loaded);
}

void gfx_set_font(gfx_font_t font)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
    gfx->plsfont(font.family, font.style, font.weight);
}

void gfx_define_colormap(const char* colormap_name,
                         int* rgba_colors,
                         size_t num_colors)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
    int_array_t* map = int_array_new();
    for (size_t i = 0; i < num_colors; ++i)
      int_array_append(map, rgba_colors[i]);
    string_ptr_unordered_map_insert_with_kv_dtors(gfx->colormaps, string_dup(colormap_name), map, string_free, DTOR(int_array_free));
  }
}

void gfx_define_palette(const char* palette_name,
                        int* rgba_colors,
                        size_t num_colors)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
    int_array_t* palette = int_array_new();
    for (size_t i = 0; i < num_colors; ++i)
      int_array_append(palette, rgba_colors[i]);
    string_ptr_unordered_map_insert_with_kv_dtors(gfx->palettes, string_dup(palette_name), palette, string_free, DTOR(int_array_free));
  }
}

void gfx_set_palette(const char* palette_name)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
    void** ptr = string_ptr_unordered_map_get(gfx->palettes, (char*)palette_name);
    if (ptr == NULL)
      polymec_error("gfx_set_palette: Unknown palette: %s", palette_name);
    int_array_t* palette = *((int_array_t**)ptr);
    int r[palette->size], g[palette->size], b[palette->size];
    double a[palette->size];
    for (size_t i = 0; i < palette->size; ++i)
    {
      r[i] = (palette->data[i] >> 24) & 255;
      g[i] = (palette->data[i] >> 16) & 255;
      b[i] = (palette->data[i] >> 8) & 255;
      a[i] = 1.0 * palette->data[i] / 255.0;
    }
    gfx->plscmap0a(r, g, b, a, (int)palette->size);
  }
}

static void gfx_figure_free(gfx_figure_t* fig)
{
  string_free(fig->title);
  string_free(fig->x_label);
  string_free(fig->y_label);
  string_free(fig->z_label);
  ptr_array_free(fig->legend);
  fig->page = NULL;
}

static gfx_figure_t* _gfx_figure_new(gfx_page_t* page)
{
  gfx_figure_t* fig = polymec_gc_malloc(sizeof(gfx_figure_t), DTOR(gfx_figure_free));
  if (page != NULL)
  {
    fig->page = page;
    fig->index = (int)page->figures->size;
  }
  else
  {
    fig->page = gfx_page_new(1, 1);
    fig->index = 1;
  }
  fig->x_min =  REAL_MAX;
  fig->x_max = -REAL_MAX;
  fig->x_scale_type = GFX_LINEAR_SCALE;
  fig->y_min =  REAL_MAX;
  fig->y_max = -REAL_MAX;
  fig->y_scale_type = GFX_LINEAR_SCALE;
  fig->title = string_dup("");
  fig->x_label = string_dup("");
  fig->y_label = string_dup("");
  fig->z_label = string_dup("");
  fig->legend = ptr_array_new();
  return fig;
}

gfx_figure_t* gfx_figure_new()
{
  return _gfx_figure_new(NULL);
}

gfx_page_t* gfx_figure_page(gfx_figure_t* fig)
{
  return fig->page;
}

void gfx_figure_set_x_label(gfx_figure_t* fig, const char* label)
{
  string_free(fig->x_label);
  fig->x_label = string_dup(label);
}

void gfx_figure_set_y_label(gfx_figure_t* fig, const char* label)
{
  string_free(fig->y_label);
  fig->y_label = string_dup(label);
}

void gfx_figure_set_z_label(gfx_figure_t* fig, const char* label)
{
  string_free(fig->z_label);
  fig->z_label = string_dup(label);
}

void gfx_figure_set_title(gfx_figure_t* fig, const char* title)
{
  string_free(fig->title);
  fig->title = string_dup(title);
}

char* gfx_figure_title(gfx_figure_t* fig)
{
  return fig->title;
}

void gfx_figure_set_x_scale(gfx_figure_t* fig, 
                            real_t x_min, real_t x_max,
                            gfx_scale_t type)
{
  ASSERT(x_min < x_max);
  fig->x_min = x_min;
  fig->x_max = x_max;
  fig->x_scale_type = type;
}

void gfx_figure_clear_x_scale(gfx_figure_t* fig,
                              gfx_scale_t type)
{
  fig->x_min = REAL_MAX;
  fig->x_max = -REAL_MAX;
  fig->x_scale_type = type;
}

void gfx_figure_set_y_scale(gfx_figure_t* fig, 
                            real_t y_min, real_t y_max,
                            gfx_scale_t type)
{
  ASSERT(y_min < y_max);
  fig->y_min = y_min;
  fig->y_max = y_max;
  fig->y_scale_type = type;
}

void gfx_figure_clear_y_scale(gfx_figure_t* fig,
                              gfx_scale_t type)
{
  fig->y_min = REAL_MAX;
  fig->y_max = -REAL_MAX;
  fig->y_scale_type = type;
}

void gfx_figure_colorbar(gfx_figure_t* fig,
                         real_t x, 
                         real_t y)
{
}

void gfx_figure_legend(gfx_figure_t* fig,
                       real_t x, 
                       real_t y)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
  }
}

static void get_data_bounds(real_t min0, real_t max0, 
                            real_t* data, size_t n,
                            double* min, double* max)
{
  *min = DBL_MAX; 
  *max = -DBL_MAX;
  if (min0 < max0) 
  {
    *min = (double)min0;
    *max = (double)max0;
  }
  else
  {
    for (size_t i = 0; i < n; ++i)
    {
      *min = MIN(min0, (double)data[i]);
      *max = MAX(max0, (double)data[i]);
    }
  }
}

static void get_xy_scale_type(gfx_figure_t* fig, int* scale_type)
{
  if (fig->x_scale_type == GFX_LINEAR_SCALE)
  {
    if (fig->y_scale_type == GFX_LINEAR_SCALE)
      *scale_type = 0;
    else // log y scaling
      *scale_type = 20;
  }
  else // log x scaling
  {
    if (fig->y_scale_type == GFX_LINEAR_SCALE)
      *scale_type = 10;
    else // log y scaling
      *scale_type = 30;
  }
}

void gfx_figure_plot(gfx_figure_t* fig, 
                     real_t* x, 
                     real_t* y, 
                     size_t n,
                     const char* glyph,
                     int color,
                     const char* label)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
    // Set the extents of the window.
    double x_min, x_max, y_min, y_max;
    get_data_bounds(fig->x_min, fig->x_max, x, n, &x_min, &x_max);
    get_data_bounds(fig->y_min, fig->y_max, y, n, &y_min, &y_max);
    int scale_type;
    get_xy_scale_type(fig, &scale_type);
    gfx->plenv0(x_min, x_max, y_min, y_max, 0, scale_type);

    // Labels and title.
    gfx->pllab(fig->x_label, fig->y_label, fig->title);

    // Set the pen color.
    gfx->plcol0(color);

    // Do the plot.
    if (glyph == NULL)
      gfx->plline((int)n, x, y);
    else
      gfx->plstring((int)n, x, y, glyph);

    // Add the label to the legend.
    legend_item_t* item = legend_item_new(label, glyph, color);
    ptr_array_append_with_dtor(fig->legend, item, DTOR(legend_item_free));
  }
}

void gfx_figure_contour(gfx_figure_t* fig)
{
}

void gfx_figure_shade(gfx_figure_t* fig)
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

static void gfx_page_free(gfx_page_t* page)
{
  gfx_t* gfx = gfx_instance();
  if (gfx->plplot != NULL)
  {
    gfx->plsstrm(page->index);
    gfx->pleop();
  }
  ptr_array_free(page->figures);
  polymec_free(page);
}

gfx_page_t* gfx_page_new(int num_rows, int num_cols)
{
  ASSERT(num_rows > 0);
  ASSERT(num_cols > 0);
  gfx_t* gfx = gfx_instance();
  gfx_page_t* page = polymec_gc_malloc(sizeof(gfx_page_t), DTOR(gfx_page_free));
  page->num_rows = num_rows;
  page->num_cols = num_cols;
  page->figures = ptr_array_new();
  for (size_t i = 0; i < num_rows * num_cols; ++i)
  {
    gfx_figure_t* fig = _gfx_figure_new(page);
    ptr_array_append(page->figures, fig);
  }
  ptr_array_append(gfx->pages, page);

  if (gfx->plplot != NULL) 
  {
    gfx->plsstrm((int)gfx->pages->size);
    gfx->plssub(num_rows, num_cols);
  }

  return page;
}

int gfx_page_num_rows(gfx_page_t* page)
{
  return page->num_rows;
}

int gfx_page_num_columns(gfx_page_t* page)
{
  return page->num_cols;
}

bool gfx_page_next(gfx_page_t* page, int* pos, gfx_figure_t** figure)
{
  return ptr_array_next(page->figures, pos, (void**)figure);
}

gfx_figure_t* gfx_page_figure(gfx_page_t* page,
                              int row,
                              int column)
{
  return page->figures->data[page->num_cols * row + column];
}

