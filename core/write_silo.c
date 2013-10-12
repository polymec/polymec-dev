// Copyright 2012-2013 Jeffrey Johnson.
// 
// This file is part of Polymec, and is licensed under the Apache License, 
// Version 2.0 (the "License"); you may not use this file except in 
// compliance with the License. You may may find the text of the license in 
// the LICENSE file at the top-level source directory, or obtain a copy of 
// it at
// 
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include "core/write_silo.h"

#undef HAVE_HDF5
#undef HAVE_TETGEN
#include "polytope_c.h"

// Create a polytope tessellation from our mesh. This tessellation 
// borrows memory from the mesh, so it shouldn't be freed.
static polytope_tessellation_t* tessellation_from_mesh(mesh_t* mesh)
{
  polytope_tessellation_t* tess = polytope_tessellation_new(3);

  // Nodes.
  tess->num_nodes = mesh->num_nodes;
  tess->nodes = malloc(sizeof(double) * 3*mesh->num_nodes);
  for (int i = 0; i < mesh->num_nodes; ++i)
  {
    tess->nodes[3*i] = mesh->nodes[i].x;
    tess->nodes[3*i+1] = mesh->nodes[i].y;
    tess->nodes[3*i+2] = mesh->nodes[i].z;
  }

  // Faces.
  tess->num_faces = mesh->num_faces;
  tess->face_offsets = malloc(sizeof(int) * (mesh->num_faces + 1));
  memcpy(tess->face_offsets, mesh->face_node_offsets, sizeof(int) * (mesh->num_faces + 1));
  tess->face_cells = malloc(sizeof(int) * 2*tess->num_faces);
  memcpy(tess->face_cells, mesh->face_cells, sizeof(int) * 2 * tess->num_faces);
  tess->face_nodes = malloc(sizeof(int) * tess->face_offsets[tess->num_faces]);
  memcpy(tess->face_nodes, mesh->face_nodes, sizeof(int) * tess->face_offsets[tess->num_faces]);

  // Cells.
  tess->num_cells = mesh->num_cells;
  tess->cell_offsets = malloc(sizeof(int) * (tess->num_cells + 1));
  memcpy(tess->cell_offsets, mesh->cell_face_offsets, sizeof(int) * (tess->num_cells + 1));
  tess->cell_faces = malloc(sizeof(int) * tess->cell_offsets[tess->num_cells]);
  memcpy(tess->cell_faces, mesh->cell_faces, sizeof(int) * tess->cell_offsets[tess->num_cells]);

  return tess;
}

void write_silo(mesh_t* mesh,
                string_ptr_unordered_map_t* node_fields,
                string_ptr_unordered_map_t* edge_fields,
                string_ptr_unordered_map_t* face_fields,
                string_ptr_unordered_map_t* cell_fields,
                const char* file_prefix,
                const char* directory,
                int cycle,
                double time,
                MPI_Comm comm,
                int num_files,
                int mpi_tag)
{
  polytope_tessellation_t* tess = tessellation_from_mesh(mesh);

  // Translate everything into polytope lingo.
  int num_node_fields = (node_fields != NULL) ? node_fields->size : 0;
  char* node_field_names[num_node_fields];
  double* node_field_data[num_node_fields];
  int pos = 0, i = 0;
  void* ptr;
  if (node_fields != NULL)
  {
    while (string_ptr_unordered_map_next(node_fields, &pos, &node_field_names[i], &ptr))
      node_field_data[i++] = ptr;
  }

  int num_edge_fields = (edge_fields != NULL) ? edge_fields->size : 0;
  char* edge_field_names[num_edge_fields];
  double* edge_field_data[num_edge_fields];
  pos = 0, i = 0;
  if (edge_fields != NULL)
  {
    while (string_ptr_unordered_map_next(edge_fields, &pos, &edge_field_names[i], &ptr))
      edge_field_data[i++] = ptr;
  }

  int num_face_fields = (face_fields != NULL) ? face_fields->size : 0;
  char* face_field_names[num_face_fields];
  double* face_field_data[num_face_fields];
  pos = 0, i = 0;
  if (face_fields != NULL)
  {
    while (string_ptr_unordered_map_next(face_fields, &pos, &face_field_names[i], &ptr))
      face_field_data[i++] = ptr;
  }

  int num_cell_fields = (cell_fields != NULL) ? cell_fields->size : 0;
  char* cell_field_names[num_cell_fields];
  double* cell_field_data[num_cell_fields];
  pos = 0, i = 0;
  if (cell_fields != NULL)
  {
    while (string_ptr_unordered_map_next(cell_fields, &pos, &cell_field_names[i], &ptr))
      cell_field_data[i++] = ptr;
  }

  polytope_write_silo(tess, num_node_fields, node_field_names, node_field_data,
                      num_edge_fields, edge_field_names, edge_field_data,
                      num_face_fields, face_field_names, face_field_data,
                      num_cell_fields, cell_field_names, cell_field_data,
                      file_prefix, directory, cycle, time,
                      comm, num_files, mpi_tag);

  polytope_tessellation_free(tess);
}

