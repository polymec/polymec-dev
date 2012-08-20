#include "geometry/create_cyl_mesh.h"
#include "geometry/block.h"
#include "geometry/create_block_mesh.h"

#ifdef __cplusplus
extern "C" {
#endif

// Mapping functions for cylindrical blocks.
typedef struct cyl_map_t
{
  double box_len;
  double axial_len;
  double radius;
  double center[3];
};

static double sqrt2 = sqrt(2.0);

static void cmapping(void* ctx, point_t* x, double* phi)
{
  cyl_map* map = (cyl_map_t*)ctx;
  phi[0] = (x[0] - 0.5)*box_len + center[0];
  phi[1] = (x[1] - 0.5)*box_len + center[1];
  phi[2] = (x[2] - 0.5)*box_len + center[2];
}

static void nmapping(void* ctx, point_t* x, double* phi)
{
  cyl_map* map = (cyl_map_t*)ctx;
  double theta = asin((x[1] - 0.5)/sqrt2);
  double r = box_len 
  phi[0] = (x[0] - 0.5)*box_len;
}

static void smapping(void* ctx, point_t* x, double* phi)
{
  cyl_map* map = (cyl_map_t*)ctx;
}

static void emapping(void* ctx, point_t* x, double* phi)
{
  cyl_map* map = (cyl_map_t*)ctx;
}

static void wmapping(void* ctx, point_t* x, double* phi)
{
  cyl_map* map = (cyl_map_t*)ctx;
}

mesh_t* create_cyl_mesh(int n_center, int n_radius, int n_axial, 
                        double box_len, double radius, double axial_len,
                        double center[])
{
  // Create 5 blocks.

  // Center block.
  int center_res[3] = {n_center, n_center, n_axial};
  cyl_map_t map = {.box_len = box_len, .axial_len = axial_len, 
                   .radius = radius, .center = center};
  sp_vtable ctable = {.eval = &cmapping};
  sp_func_t* cmap = sp_func_new("center", map, ctable, SP_INHOMOGENEOUS, 3);
  block_t* c = block_new(center_res, cmap);

  // North block.
  int north_res[3] = {n_center, n_radius, n_axial};
  sp_vtable ntable = {.eval = &nmapping};
  sp_func_t* nmap = sp_func_new("north", map, ntable, SP_INHOMOGENEOUS, 3);
  block_t* n = block_new(north_res, nmap);

  // South block.
  int south_res[3] = {n_center, n_radius, n_axial};
  sp_vtable stable = {.eval = &smapping};
  sp_func_t* smap = sp_func_new("south", map, stable, SP_INHOMOGENEOUS, 3);
  block_t* s = block_new(south_res, smap);

  // East block.
  int east_res[3] = {n_radius, n_center, n_axial};
  sp_vtable etable = {.eval = &emapping};
  sp_func_t* emap = sp_func_new("east", map, etable, SP_INHOMOGENEOUS, 3);
  block_t* e = block_new(east_res, emap);

  // West block.
  int west_res[3] = {n_radius, n_center, n_axial};
  sp_vtable wtable = {.eval = &wmapping};
  sp_func_t* wmap = sp_func_new("west", map, wtable, SP_INHOMOGENEOUS, 3);
  block_t* w = block_new(west_res, wmap);

  // Attach the blocks to one another.
  block_attach(c, 3, n, 2);
  block_attach(c, 2, s, 3);
  block_attach(c, 1, e, 0);
  block_attach(c, 0, w, 1);

  block_attach(n, 0, w, 3);
  block_attach(w, 2, s, 0);
  block_attach(s, 1, e, 2);
  block_attach(e, 3, n, 1);

  // Use the blocks to assemble a mesh.
  block_t* blocks[5] = {c, n, s, e, w};
  mesh_t* mesh = create_block_mesh(5, blocks);

  // Clean up.
  block_free(n);
  block_free(s);
  block_free(e);
  block_free(w);

  return mesh;
}

#ifdef __cplusplus
}
#endif

