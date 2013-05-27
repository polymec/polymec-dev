#include <stdlib.h>
#include <string.h>
#include "core/edit_mesh.h"
#include "core/mesh_storage.h"

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
    // We have exceeded our current node capacity. So we have to 
    // reallocate memory.
    while (mesh->num_nodes+1 > mesh->storage->node_capacity)
      mesh->storage->node_capacity *= 2;

    // Because the addresses of the nodes within the mesh are invalidated 
    // by a realloc, we allocate a new chunk of memory for the nodes and 
    // copy the data from the old chunk to the new, and then refresh the 
    // other mesh elements w.r.t. these new pointers.
    node_t* new_nodes = ARENA_MALLOC(mesh->arena, sizeof(node_t)*mesh->storage->node_capacity, 0);
    memcpy(new_nodes, mesh->nodes, mesh->num_nodes * sizeof(node_t));
    memset(new_nodes + mesh->num_nodes, 0, sizeof(node_t)*(mesh->storage->node_capacity - mesh->num_nodes));

    // Refresh each edge's node pointers.
    for (int e = 0; e < mesh->num_edges; ++e)
    {
      edge_t* edge = &mesh->edges[e];
      if (edge->node1 != NULL)
      {
        ptrdiff_t offset = edge->node1 - &mesh->nodes[0];
        edge->node1 = new_nodes + offset;
      }
      if (edge->node2 != NULL)
      {
        ptrdiff_t offset = edge->node2 - &mesh->nodes[0];
        edge->node2 = new_nodes + offset;
      }
    }

    // Get rid of the old nodes.
    ARENA_FREE(mesh->arena, mesh->nodes);
    mesh->nodes = new_nodes;
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
    mesh->num_nodes -= 1;
  }
}

int mesh_add_edge(mesh_t* mesh)
{
  if (mesh->num_edges+1 > mesh->storage->edge_capacity)
  {
    while (mesh->num_edges+1 > mesh->storage->edge_capacity)
      mesh->storage->edge_capacity *= 2;

    // Because the addresses of the edges within the mesh are invalidated 
    // by a realloc, we allocate a new chunk of memory for the edges and 
    // copy the data from the old chunk to the new, and then refresh the 
    // other mesh elements w.r.t. these new pointers.
    edge_t* new_edges = ARENA_MALLOC(mesh->arena, sizeof(edge_t)*mesh->storage->edge_capacity, 0);
    memcpy(new_edges, mesh->edges, mesh->num_edges * sizeof(edge_t));
    memset(new_edges + mesh->num_edges, 0, sizeof(edge_t)*(mesh->storage->edge_capacity-mesh->num_edges));

    // Refresh each face's edge pointers.
    for (int f = 0; f < mesh->num_faces; ++f)
    {
      face_t* face = &mesh->faces[f];
      for (int e = 0; e < face->num_edges; ++e)
      {
        ptrdiff_t offset = face->edges[e] - &mesh->edges[0];
        face->edges[e] = new_edges + offset;
      }
    }

    // Get rid of the old edges.
    ARENA_FREE(mesh->arena, mesh->edges);
    mesh->edges = new_edges;
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
    mesh->num_edges -= 1;
  }
}

int mesh_add_face(mesh_t* mesh)
{
  if (mesh->num_faces+1 > mesh->storage->face_capacity)
  {
    while (mesh->num_faces+1 > mesh->storage->face_capacity)
      mesh->storage->face_capacity *= 2;

    // Because the addresses of the faces within the mesh are invalidated 
    // by a realloc, we allocate a new chunk of memory for the faces and 
    // copy the data from the old chunk to the new, and then refresh the 
    // other mesh elements w.r.t. these new pointers.
    face_t* new_faces = ARENA_MALLOC(mesh->arena, sizeof(face_t)*mesh->storage->face_capacity, 0);
    memcpy(new_faces, mesh->faces, mesh->num_faces * sizeof(face_t));
    memset(new_faces + mesh->num_faces, 0, sizeof(face_t)*(mesh->storage->face_capacity-mesh->num_faces));

    // Refresh each cell's face pointers.
    for (int c = 0; c < mesh->num_cells; ++c)
    {
      cell_t* cell = &mesh->cells[c];
      for (int f = 0; f < cell->num_faces; ++f)
      {
        ptrdiff_t offset = cell->faces[f] - &mesh->faces[0];
        cell->faces[f] = new_faces + offset;
      }
    }

    // Get rid of the old faces.
    ARENA_FREE(mesh->arena, mesh->faces);
    mesh->faces = new_faces;
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
    mesh->num_faces -= 1;
  }
}

int mesh_add_cell(mesh_t* mesh)
{
  if (mesh->num_cells+1 > mesh->storage->cell_capacity)
  {
    while (mesh->num_cells+1 > mesh->storage->cell_capacity)
      mesh->storage->cell_capacity *= 2;

    // Because the addresses of the cells within the mesh are invalidated 
    // by a realloc, we allocate a new chunk of memory for the cells and 
    // copy the data from the old chunk to the new, and then refresh the 
    // other mesh elements w.r.t. these new pointers.
    cell_t* new_cells = ARENA_MALLOC(mesh->arena, sizeof(cell_t)*mesh->storage->cell_capacity, 0);
    memcpy(new_cells, mesh->cells, mesh->num_cells * sizeof(cell_t));
    memset(new_cells + mesh->num_cells, 0, sizeof(cell_t)*(mesh->storage->cell_capacity-mesh->num_cells));

    // Refresh each face's cell pointers.
    for (int f = 0; f < mesh->num_faces; ++f)
    {
      face_t* face = &mesh->faces[f];
      if (face->cell1 != NULL)
      {
        ptrdiff_t offset = face->cell1 - &mesh->cells[0];
        face->cell1 = new_cells + offset;
      }
      if (face->cell2 != NULL)
      {
        ptrdiff_t offset = face->cell2 - &mesh->cells[0];
        face->cell2 = new_cells + offset;
      }
    }

    // Get rid of the old cells.
    ARENA_FREE(mesh->arena, mesh->cells);
    mesh->cells = new_cells;
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
    mesh->num_cells -= 1;
  }
}

void mesh_attach_edge_to_face(mesh_t* mesh, edge_t* edge, face_t* face)
{
  ASSERT((edge - &mesh->edges[0]) >= 0);
  ASSERT((edge - &mesh->edges[0]) < mesh->num_edges);
  ASSERT((face - &mesh->faces[0]) >= 0);
  ASSERT((face - &mesh->faces[0]) < mesh->num_faces);

  if (face->edges == NULL)
  {
    face->num_edges = 1;
    face->edges = ARENA_MALLOC(mesh->arena, sizeof(edge_t*)*4, 0);
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
    face->edges = ARENA_REALLOC(mesh->arena, face->edges, sizeof(edge_t*)*ne, 0);
    face->edges[face->num_edges-1] = edge;
  }
}

void mesh_attach_face_to_cell(mesh_t* mesh, face_t* face, cell_t* cell)
{
  ASSERT((face - &mesh->faces[0]) >= 0);
  ASSERT((face - &mesh->faces[0]) < mesh->num_faces);
  ASSERT((cell - &mesh->cells[0]) >= 0);
  ASSERT((cell - &mesh->cells[0]) < mesh->num_cells);

  // Make sure this face isn't already attached to two cells.
  ASSERT((face->cell1 == NULL) || (face->cell2 == NULL));

  if (cell->faces == NULL)
  {
    cell->num_faces = 1;
    cell->faces = ARENA_MALLOC(mesh->arena, sizeof(face_t*)*4, 0);
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
    cell->faces = ARENA_REALLOC(mesh->arena, cell->faces, sizeof(face_t*)*nf, 0);
    cell->faces[cell->num_faces-1] = face;
  }
  if (face->cell1 == NULL)
    face->cell1 = cell;
  else
    face->cell2 = cell;
}

void mesh_detach_edge_from_face(mesh_t* mesh, edge_t* edge, face_t* face)
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

void mesh_detach_face_from_cell(mesh_t* mesh, face_t* face, cell_t* cell)
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

