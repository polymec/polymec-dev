// Copyright (c) 2012-2017, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_GFX_H
#define POLYMEC_GFX_H

#include "core/polymec.h"

// This class represents an interactive 2D figure. Figures are garbage 
// collected.
typedef struct gfx_figure_t gfx_figure_t;

// This class represents a page containing one or more figures. Pages are 
// garbage collected.
typedef struct gfx_page_t gfx_page_t;

// Font families.
typedef enum
{
  GFX_FONT_SANS,
  GFX_FONT_SERIF,
  GFX_FONT_MONO,
  GFX_FONT_SCRIPT,
  GFX_FONT_SYMBOL
} gfx_font_family_t;

// Font styles.
typedef enum
{
  GFX_FONT_UPRIGHT,
  GFX_FONT_ITALIC,
  GFX_FONT_OBLIQUE
} gfx_font_style_t;

// Font weights.
typedef enum
{
  GFX_FONT_MEDIUM,
  GFX_FONT_BOLD
} gfx_font_weight_t;

// Creates a new figure on its own page.
gfx_figure_t* gfx_figure_new();

// Returns an internal pointer to the page for the given figure.
gfx_page_t* gfx_figure_page(gfx_figure_t* fig);

// Sets the label for the x axis of the given figure.
void gfx_figure_set_x_label(gfx_figure_t* fig, const char* label);

// Sets the label for the y axis of the given figure.
void gfx_figure_set_y_label(gfx_figure_t* fig, const char* label);

// Sets the label for the z axis of the given figure, if it has one. 
// If the figure has no z label, this function has no effect.
void gfx_figure_set_z_label(gfx_figure_t* fig, const char* label);

// Sets the title of the figure.
void gfx_figure_set_title(gfx_figure_t* fig, const char* title);

// Sets the discrete color map for the figure by specifying a given 
// set of RGBA colors.
void gfx_figure_set_discrete_colors(gfx_figure_t* fig,
                                    int* rgba,
                                    size_t num_colors);

// Sets the continuous color map for the figure by specifying a given 
// set of RGBA colors.
void gfx_figure_set_cont_colors(gfx_figure_t* fig,
                                int* rgba,
                                size_t num_points);

// Sets the font to use for subsequent text, specifying the family, style, 
// and weight.
void gfx_figure_set_font(gfx_figure_t* fig,
                         gfx_font_family_t family,
                         gfx_font_style_t style,
                         gfx_font_weight_t weight);

// Adds a color bar to the figure, or resets the existing one.
void gfx_figure_colorbar(gfx_figure_t* fig);

// Adds a legend to the figure, or resets the existing one.
void gfx_figure_legend(gfx_figure_t* fig);

// Generates a curve consisting of n (x, y) points into the given figure. If 
// label is non-NULL, the plot will be annotated in any associated legend.
void gfx_figure_plot(gfx_figure_t* fig, 
                     real_t* x, 
                     real_t* y, 
                     size_t n,
                     char* label);

// Generates a scatter plot of n (x, y) point data in the given figure. 
void gfx_figure_scatter(gfx_figure_t* fig,
                        real_t* x, 
                        real_t* y,
                        size_t n);

// Generates a contour plot in the given figure.
void gfx_figure_contour(gfx_figure_t* fig);

// Generates a surface plot in the given figure.
void gfx_figure_surface(gfx_figure_t* fig);

// Generates a vector plot in the given figure.
void gfx_figure_quiver(gfx_figure_t* fig);

// Generates an image in the given figure.
void gfx_figure_image(gfx_figure_t* fig);

// Clears the given figure.
void gfx_figure_clear(gfx_figure_t* fig);

// Creates a new page that contains figures in each of its rows and columns.
gfx_page_t* gfx_page_new(int num_rows, int num_cols);

// Returns the number of rows in the given page.
int gfx_page_num_rows(gfx_page_t* page);

// Returns the number of columns in the given page.
int gfx_page_num_columns(gfx_page_t* page);

// Traverses the figures in the given page. Set pos to zero to reset the 
// traversal.
bool gfx_page_next(gfx_page_t* page, int* pos, gfx_figure_t** figure);

// Returns the figure in the given row and column within the given page.
gfx_figure_t* gfx_page_figure(gfx_page_t* page,
                              int row,
                              int column);

#endif
