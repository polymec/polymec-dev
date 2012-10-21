#include <stdlib.h>
#include <string.h>
#include "core/edit_mesh.h"
#include "core/mesh_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

#define MAX(x, y) ((x >= y) ? x : y)

// This function rounds the given number up to the nearest power of 2.
static int round_to_pow2(int x)
{
  int y = 2;
  while (y < x) y *= 2;
  return y;
}

int mesh_add_node(mesh_t* mesh)
{
  if (mesh->num_nodes+1 > mesh->storage->node_capacity)
  {
    while (mesh->num_nodes+1 > mesh->storage->node_capacity)
      mesh->storage->node_capacity *= 2;
    mesh->nodes = arena_realloc(mesh->arena, mesh->nodes, sizeof(node_t)*mesh->storage->node_capacity, 0);
  }
  mesh->num_nodes++;
  return mesh->num_nodes-1;
}

void mesh_delete_node(mesh_t* mesh, int i)
{
  // Swap the ith node with the end.
  if (i < mesh->num_nodes)
  {
    int last = mesh->num_nodes - 1;
    node_t tmp = mesh->nodes[last];
    mesh->nodes[last] = mesh->nodes[i];
    mesh->nodes[i] = tmp;
  }
}

int mesh_add_edge(mesh_t* mesh)
{
  if (mesh->num_edges+1 > mesh->storage->edge_capacity)
  {
    int old_cap = mesh->storage->edge_capacity;
    while (mesh->num_edges+1 > mesh->storage->edge_capacity)
      mesh->storage->edge_capacity *= 2;
    mesh->edges = arena_realloc(mesh->arena, mesh->edges, sizeof(edge_t)*mesh->storage->edge_capacity, 0);
    memset(mesh->edges + old_cap, 0, sizeof(edge_t)*(mesh->storage->edge_capacity-old_cap));
  }
  mesh->num_edges++;
  return mesh->num_edges-1;
}

void mesh_delete_edge(mesh_t* mesh, int i)
{
  // Swap the ith edge with the end.
  if (i < mesh->num_edges)
  {
    int last = mesh->num_edges - 1;
    edge_t tmp = mesh->edges[last];
    mesh->edges[last] = mesh->edges[i];
    mesh->edges[i] = tmp;
  }
}

int mesh_add_face(mesh_t* mesh)
{
  if (mesh->num_faces+1 > mesh->storage->face_capacity)
  {
    int old_cap = mesh->storage->face_capacity;
    while (mesh->num_faces+1 > mesh->storage->face_capacity)
      mesh->storage->face_capacity *= 2;
    mesh->faces = arena_realloc(mesh->arena, mesh->faces, sizeof(face_t)*mesh->storage->face_capacity, 0);
    memset(mesh->faces + old_cap, 0, sizeof(face_t)*(mesh->storage->face_capacity-old_cap));
  }
  mesh->num_faces++;
  return mesh->num_faces-1;
}

void mesh_delete_face(mesh_t* mesh, int i)
{
  // Swap the ith face with the end.
  if (i < mesh->num_faces)
  {
    int last = mesh->num_faces - 1;
    face_t tmp = mesh->faces[last];
    mesh->faces[last] = mesh->faces[i];
    mesh->faces[i] = tmp;
  }
}

int mesh_add_cell(mesh_t* mesh)
{
  if (mesh->num_cells+1 > mesh->storage->cell_capacity)
  {
    int old_cap = mesh->storage->cell_capacity;
    while (mesh->num_cells+1 > mesh->storage->cell_capacity)
      mesh->storage->cell_capacity *= 2;
    mesh->cells = arena_realloc(mesh->arena, mesh->cells, sizeof(cell_t)*mesh->storage->cell_capacity, 0);
    memset(mesh->cells + old_cap, 0, sizeof(cell_t)*(mesh->storage->cell_capacity-old_cap));
  }
  mesh->num_cells++;
  return mesh->num_cells-1;
}

void mesh_delete_cell(mesh_t* mesh, int i)
{
  // Swap the ith cell with the end.
  if (i < mesh->num_cells)
  {
    int last = mesh->num_cells - 1;
    cell_t tmp = mesh->cells[last];
    mesh->cells[last] = mesh->cells[i];
    mesh->cells[i] = tmp;
  }
}

void mesh_add_edge_to_face(mesh_t* mesh, edge_t* edge, face_t* face)
{
  if (face->edges == NULL)
  {
    face->num_edges = 1;
    face->edges = arena_malloc(mesh->arena, sizeof(edge_t*)*4, 0);
    face->edges[0] = edge;
  }
  else
  {
#ifndef NDEBUG
    // Make sure this edge isn't already attached.
    for (int e = 0; e < face->num_edges; ++e)
    {
      ASSERT(edge != face->edges[e]);
    }
#endif
    face->num_edges++;
    int ne = MAX(round_to_pow2(face->num_edges), 4);
    face->edges = arena_realloc(mesh->arena, face->edges, sizeof(edge_t*)*ne, 0);
    face->edges[face->num_edges-1] = edge;
  }
}

void mesh_add_face_to_cell(mesh_t* mesh, face_t* face, cell_t* cell)
{
  // Make sure this face isn't already attached to two cells.
  ASSERT((face->cell1 == NULL) || (face->cell2 == NULL));

  if (cell->faces == NULL)
  {
    cell->num_faces = 1;
    cell->faces = arena_malloc(mesh->arena, sizeof(face_t*)*4, 0);
    cell->faces[0] = face;
  }
  else
  {
#ifndef NDEBUG
    // Make sure this face isn't already attached.
    for (int f = 0; f < cell->num_faces; ++f)
    {
      ASSERT(face != cell->faces[f]);
    }
#endif
    cell->num_faces++;
    int nf = MAX(round_to_pow2(cell->num_faces), 4);
    cell->faces = arena_realloc(mesh->arena, cell->faces, sizeof(face_t*)*nf, 0);
    cell->faces[cell->num_faces-1] = face;
  }
  if (face->cell1 == NULL)
    face->cell1 = cell;
  else
    face->cell2 = cell;
}

void mesh_remove_edge_from_face(mesh_t* mesh, edge_t* edge, face_t* face)
{
  for (int e = 0; e < face->num_edges; ++e)
  {
    if (face->edges[e] == edge)
    {
      face->num_edges--;
      face->edges[e] = face->edges[face->num_edges];
      break;
    }
  }
}

void mesh_remove_face_from_cell(mesh_t* mesh, face_t* face, cell_t* cell)
{
  for (int f = 0; f < cell->num_faces; ++f)
  {
    if (cell->faces[f] == face)
    {
      cell->num_faces--;
      cell->faces[f] = cell->faces[cell->num_faces];
      break;
    }
  }
}

#ifdef __cplusplus
}
#endif

