#include <gc/gc.h>
#include "core/mesh_delta.h"
#include "core/edit_mesh.h"

// Function prototypes for the base class.
typedef void (*mesh_delta_apply_func)(void*, mesh_t*);
typedef mesh_delta_t* (*mesh_delta_inverse_func)(void*);
typedef void (*mesh_delta_dtor)(void*);

typedef struct 
{
  mesh_delta_apply_func apply;
  mesh_delta_inverse_func inverse;
  mesh_delta_dtor dtor;
} mesh_delta_vtable;

// Base class definition.
struct mesh_delta_t 
{
  char* name;
  void* context;
  mesh_delta_vtable vtable;
};

static void mesh_delta_free(void* ctx, void* dummy)
{
  mesh_delta_t* delta = ctx;
  free(delta->name);
  if ((delta->context != NULL) && (delta->vtable.dtor))
    delta->vtable.dtor(delta->context);
}

// Base class constructor.
static mesh_delta_t* mesh_delta_new(const char* name, 
                                    void* context,
                                    mesh_delta_vtable vtable)
{
  ASSERT(name != NULL);
  mesh_delta_t* delta = GC_MALLOC(sizeof(mesh_delta_t));
  delta->name = strdup(name);
  delta->context = context;
  delta->vtable = vtable;
  GC_register_finalizer(delta, mesh_delta_free, delta, NULL, NULL);
  return delta;
}

void mesh_delta_apply(mesh_delta_t* delta, mesh_t* mesh)
{
  delta->vtable.apply(delta->context, mesh);
}

mesh_delta_t* mesh_delta_inverse(mesh_delta_t* delta)
{
  return delta->vtable.inverse(delta->context);
}

void mesh_delta_fprintf(mesh_delta_t* delta, FILE* file)
{
  fprintf(file, "%s", delta->name);
}

// ---------------------------- Swap delta ------------------------------
typedef struct 
{
  mesh_centering_t type;
  void (*swap)(mesh_t* mesh, int, int);
  int index1, index2;
} swap_delta_t;

static void swap_apply(void* context, mesh_t* mesh)
{
  swap_delta_t* swap = context;
  swap->swap(mesh, swap->index1, swap->index2);
}

static mesh_delta_t* swap_inverse(void* context)
{
  swap_delta_t* swap = context;
  return swap_mesh_delta_new(swap->type, swap->index1, swap->index2);
}

static void swap_node(mesh_t* mesh, int index1, int index2)
{
  ASSERT(index1 < mesh->num_nodes);
  ASSERT(index2 < mesh->num_nodes);
  node_t tmp = mesh->nodes[index2];
  mesh->nodes[index2] = mesh->nodes[index1];
  mesh->nodes[index1] = tmp;
}

static void swap_edge(mesh_t* mesh, int index1, int index2)
{
  ASSERT(index1 < mesh->num_edges);
  ASSERT(index2 < mesh->num_edges);
  edge_t tmp = mesh->edges[index2];
  mesh->edges[index2] = mesh->edges[index1];
  mesh->edges[index1] = tmp;
}

static void swap_face(mesh_t* mesh, int index1, int index2)
{
  ASSERT(index1 < mesh->num_faces);
  ASSERT(index2 < mesh->num_faces);
  face_t tmp = mesh->faces[index2];
  mesh->faces[index2] = mesh->faces[index1];
  mesh->faces[index1] = tmp;
}

static void swap_cell(mesh_t* mesh, int index1, int index2)
{
  ASSERT(index1 < mesh->num_cells);
  ASSERT(index2 < mesh->num_cells);
  cell_t tmp = mesh->cells[index2];
  mesh->cells[index2] = mesh->cells[index1];
  mesh->cells[index1] = tmp;
}

mesh_delta_t* swap_mesh_delta_new(mesh_centering_t type, int index1, int index2)
{
  ASSERT(index1 >= 0);
  ASSERT(index2 >= 0);
  ASSERT(index1 != index2);
  swap_delta_t* swap = malloc(sizeof(swap_delta_t));
  swap->type = type;
  char name[1024];
  switch (type)
  {
    case MESH_NODE:
      snprintf(name, 1024, "swap nodes %d and %d", index1, index2);
      swap->swap = swap_node;
      break;
    case MESH_EDGE:
      snprintf(name, 1024, "swap edges %d and %d", index1, index2);
      swap->swap = swap_edge;
      break;
    case MESH_FACE:
      snprintf(name, 1024, "swap faces %d and %d", index1, index2);
      swap->swap = swap_face;
      break;
    case MESH_CELL:
      snprintf(name, 1024, "swap cells %d and %d", index1, index2);
      swap->swap = swap_cell;
      break;
  }
  swap->index1 = index1;
  swap->index2 = index2;
  mesh_delta_vtable vtable = {.apply = swap_apply, .inverse = swap_inverse};
  return mesh_delta_new(name, swap, vtable);
}

// ---------------------------- Append delta ------------------------------
typedef struct 
{
  mesh_centering_t type;
  void (*append)(mesh_t* mesh, point_t* coord);
  point_t coord;
} append_delta_t;

static void append_apply(void* context, mesh_t* mesh)
{
  append_delta_t* append = context;
  append->append(mesh, &append->coord);
}

static mesh_delta_t* append_inverse(void* context)
{
  append_delta_t* append = context;
  return pop_mesh_delta_new(append->type);
}

static void append_node(mesh_t* mesh, point_t* coord)
{
  mesh_add_node(mesh);
  mesh->nodes[mesh->num_nodes-1].x = coord->x;
  mesh->nodes[mesh->num_nodes-1].y = coord->y;
  mesh->nodes[mesh->num_nodes-1].z = coord->z;
}

static void append_edge(mesh_t* mesh, point_t* coord)
{
  mesh_add_edge(mesh);
}

static void append_face(mesh_t* mesh, point_t* coord)
{
  mesh_add_face(mesh);
}

static void append_cell(mesh_t* mesh, point_t* coord)
{
  mesh_add_cell(mesh);
}

mesh_delta_t* append_mesh_delta_new(mesh_centering_t type)
{
  ASSERT(type != MESH_NODE); // This is handled elsewhere.
  append_delta_t* append = malloc(sizeof(append_delta_t));
  append->type = type;
  char name[1024];
  switch (type)
  {
    case MESH_EDGE:
      snprintf(name, 1024, "append edge to mesh");
      append->append = append_edge;
      break;
    case MESH_FACE:
      snprintf(name, 1024, "append face to mesh");
      append->append = append_face;
      break;
    case MESH_CELL:
      snprintf(name, 1024, "append cell to mesh");
      append->append = append_cell;
      break;
    default:
      ASSERT(false);
      break;
  }
  mesh_delta_vtable vtable = {.apply = append_apply, .inverse = append_inverse};
  return mesh_delta_new(name, append, vtable);
}

mesh_delta_t* append_node_mesh_delta_new(point_t* x)
{
  append_delta_t* append = malloc(sizeof(append_delta_t));
  append->append = append_node;
  append->type = MESH_NODE;
  append->coord.x = x->x;
  append->coord.y = x->y;
  append->coord.z = x->z;
  char name[1024];
  snprintf(name, 1024, "append node to mesh");
  mesh_delta_vtable vtable = {.apply = append_apply, .inverse = append_inverse};
  return mesh_delta_new(name, append, vtable);
}

// ---------------------------- Pop delta ------------------------------
typedef struct 
{
  mesh_centering_t type;
  void (*pop)(mesh_t* mesh);
} pop_delta_t;

static void pop_apply(void* context, mesh_t* mesh)
{
  pop_delta_t* pop = context;
  pop->pop(mesh);
}

static mesh_delta_t* pop_inverse(void* context)
{
  pop_delta_t* pop = context;
  return append_mesh_delta_new(pop->type);
}

static void pop_node(mesh_t* mesh)
{
  mesh_delete_node(mesh, mesh->num_nodes-1);
}

static void pop_edge(mesh_t* mesh)
{
  mesh_delete_edge(mesh, mesh->num_edges-1);
}

static void pop_face(mesh_t* mesh)
{
  mesh_delete_face(mesh, mesh->num_faces-1);
}

static void pop_cell(mesh_t* mesh)
{
  mesh_delete_cell(mesh, mesh->num_cells-1);
}

mesh_delta_t* pop_mesh_delta_new(mesh_centering_t type)
{
  pop_delta_t* pop = malloc(sizeof(pop_delta_t));
  pop->type = type;
  char name[1024];
  switch (type)
  {
    case MESH_NODE:
      snprintf(name, 1024, "pop node from mesh");
      pop->pop = pop_node;
      break;
    case MESH_EDGE:
      snprintf(name, 1024, "pop edge from mesh");
      pop->pop = pop_edge;
      break;
    case MESH_FACE:
      snprintf(name, 1024, "pop face from mesh");
      pop->pop = pop_face;
      break;
    case MESH_CELL:
      snprintf(name, 1024, "pop cell from mesh");
      pop->pop = pop_cell;
      break;
  }
  mesh_delta_vtable vtable = {.apply = pop_apply, .inverse = pop_inverse};
  return mesh_delta_new(name, pop, vtable);
}

// ---------------------------- Attach delta ------------------------------
typedef struct 
{
  mesh_centering_t type;
  void (*attach)(mesh_t* mesh, int, int);
  int index, parent_index;
} attach_delta_t;

static void attach_apply(void* context, mesh_t* mesh)
{
  attach_delta_t* attach = context;
  attach->attach(mesh, attach->index, attach->parent_index);
}

static mesh_delta_t* attach_inverse(void* context)
{
  attach_delta_t* attach = context;
  return detach_mesh_delta_new(attach->type, attach->index, attach->parent_index);
}

static void attach_node(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_nodes);
  ASSERT(parent_index < mesh->num_edges);
  edge_t* parent = &mesh->edges[parent_index];
  node_t* node = &mesh->nodes[index];
  if (parent->node1 == NULL)
    parent->node1 = node;
  else if (parent->node2 == NULL)
    parent->node2 = node;
}

static void attach_edge(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_edges);
  ASSERT(parent_index < mesh->num_faces);
  mesh_attach_edge_to_face(mesh, &mesh->edges[index], &mesh->faces[parent_index]);
}

static void attach_face(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_faces);
  ASSERT(parent_index < mesh->num_cells);
  mesh_attach_face_to_cell(mesh, &mesh->faces[index], &mesh->cells[parent_index]);
}

mesh_delta_t* attach_mesh_delta_new(mesh_centering_t type, int index, int parent_index)
{
  ASSERT(type != MESH_CELL);
  ASSERT(index >= 0);
  ASSERT(parent_index >= 0);
  attach_delta_t* attach = malloc(sizeof(attach_delta_t));
  attach->type = type;
  char name[1024];
  switch (type)
  {
    case MESH_NODE:
      snprintf(name, 1024, "attach node %d to edge %d", index, parent_index);
      attach->attach = attach_node;
      break;
    case MESH_EDGE:
      snprintf(name, 1024, "attach edge %d to face %d", index, parent_index);
      attach->attach = attach_edge;
      break;
    case MESH_FACE:
      snprintf(name, 1024, "attach face %d to cell %d", index, parent_index);
      attach->attach = attach_face;
      break;
    default:
      // Shouldn't get here!
      break;
  }
  attach->index = index;
  attach->parent_index = parent_index;
  mesh_delta_vtable vtable = {.apply = attach_apply, .inverse = attach_inverse};
  return mesh_delta_new(name, attach, vtable);
}

// ---------------------------- Detach delta ------------------------------
typedef struct 
{
  mesh_centering_t type;
  void (*detach)(mesh_t* mesh, int, int);
  int index, parent_index;
} detach_delta_t;

static void detach_apply(void* context, mesh_t* mesh)
{
  detach_delta_t* detach = context;
  detach->detach(mesh, detach->index, detach->parent_index);
}

static mesh_delta_t* detach_inverse(void* context)
{
  detach_delta_t* detach = context;
  return attach_mesh_delta_new(detach->type, detach->index, detach->parent_index);
}

static void detach_node(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_nodes);
  ASSERT(parent_index < mesh->num_edges);
  node_t* node = &mesh->nodes[index];
  edge_t* parent = &mesh->edges[parent_index];
  if (parent->node1 == node)
    parent->node1 = NULL;
  else if (parent->node2 == node)
    parent->node2 = NULL;
}

static void detach_edge(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_edges);
  ASSERT(parent_index < mesh->num_faces);
  mesh_detach_edge_from_face(mesh, &mesh->edges[index], &mesh->faces[parent_index]);
}

static void detach_face(mesh_t* mesh, int index, int parent_index)
{
  ASSERT(index < mesh->num_faces);
  ASSERT(parent_index < mesh->num_cells);
  mesh_detach_face_from_cell(mesh, &mesh->faces[index], &mesh->cells[parent_index]);
}

mesh_delta_t* detach_mesh_delta_new(mesh_centering_t type, int index, int parent_index)
{
  ASSERT(type != MESH_CELL);
  ASSERT(index >= 0);
  ASSERT(parent_index >= 0);
  detach_delta_t* detach = malloc(sizeof(detach_delta_t));
  detach->type = type;
  char name[1024];
  switch (type)
  {
    case MESH_NODE:
      snprintf(name, 1024, "detach node %d from edge %d", index, parent_index);
      detach->detach = detach_node;
      break;
    case MESH_EDGE:
      snprintf(name, 1024, "detach edge %d from face %d", index, parent_index);
      detach->detach = detach_edge;
      break;
    case MESH_FACE:
      snprintf(name, 1024, "detach face %d from cell %d", index, parent_index);
      detach->detach = detach_face;
      break;
    default:
      // Shouldn't get here!
      break;
  }
  detach->index = index;
  detach->parent_index = parent_index;
  mesh_delta_vtable vtable = {.apply = detach_apply, .inverse = detach_inverse};
  return mesh_delta_new(name, detach, vtable);
}

