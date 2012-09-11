#include <stdlib.h>
#include "core/edit_mesh.h"
#include "core/mesh_storage.h"

#ifdef __cplusplus
extern "C" {
#endif

int add_mesh_node(mesh_t* mesh)
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

void delete_mesh_node(mesh_t* mesh, int i)
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

int add_mesh_edge(mesh_t* mesh)
{
  if (mesh->num_edges+1 > mesh->storage->edge_capacity)
  {
    while (mesh->num_edges+1 > mesh->storage->edge_capacity)
      mesh->storage->edge_capacity *= 2;
    mesh->edges = arena_realloc(mesh->arena, mesh->edges, sizeof(edge_t)*mesh->storage->edge_capacity, 0);
  }
  mesh->num_edges++;
  return mesh->num_edges-1;
}

void delete_mesh_edge(mesh_t* mesh, int i)
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

int add_mesh_face(mesh_t* mesh)
{
  if (mesh->num_faces+1 > mesh->storage->face_capacity)
  {
    while (mesh->num_faces+1 > mesh->storage->face_capacity)
      mesh->storage->face_capacity *= 2;
    mesh->faces = arena_realloc(mesh->arena, mesh->faces, sizeof(face_t)*mesh->storage->face_capacity, 0);
  }
  mesh->num_faces++;
  return mesh->num_faces-1;
}

void delete_mesh_face(mesh_t* mesh, int i)
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

int add_mesh_cell(mesh_t* mesh)
{
  if (mesh->num_cells+1 > mesh->storage->cell_capacity)
  {
    while (mesh->num_cells+1 > mesh->storage->cell_capacity)
      mesh->storage->cell_capacity *= 2;
    mesh->cells = arena_realloc(mesh->arena, mesh->cells, sizeof(cell_t)*mesh->storage->cell_capacity, 0);
  }
  mesh->num_cells++;
  return mesh->num_cells-1;
}

void delete_mesh_cell(mesh_t* mesh, int i)
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

void add_face_edge(face_t* face, edge_t* edge)
{
}

void add_cell_face(cell_t* cell, face_t* face)
{
}

void remove_face_edge(face_t* face, edge_t* edge)
{
}

void remove_cell_face(cell_t* cell, face_t* face)
{
}

#ifdef __cplusplus
}
#endif

