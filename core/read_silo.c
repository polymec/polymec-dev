// Copyright (c) 2012-2014, Jeffrey N. Johnson
// All rights reserved.
// 
// Redistribution and use in source and binary forms, with or without 
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this 
// list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice, 
// this list of conditions and the following disclaimer in the documentation 
// and/or other materials provided with the distribution.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE 
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE 
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE 
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "core/read_silo.h"
#include "polytope_c.h"

// Create a mesh from a polytope tessellation. 
static mesh_t* mesh_from_tessellation(polytope_tessellation_t* tess)
{
  // FIXME
  return NULL;
}

static void field_map_kv_dtor(char* key, void* value)
{
  free(key);
  free(value);
}

void read_silo_mesh(mesh_t** mesh,
                    string_ptr_unordered_map_t** fields,
                    const char* file_prefix,
                    const char* directory,
                    int cycle,
                    real_t* time,
                    MPI_Comm comm,
                    int num_files,
                    int mpi_tag)
{
  // Read the stuff from the file.
  int num_fields, num_node_tags, num_edge_tags, num_face_tags, num_cell_tags; 
  int *node_tag_sizes, *edge_tag_sizes, *face_tag_sizes, *cell_tag_sizes; 
  int **node_tags, **edge_tags, **face_tags, **cell_tags;
  char **field_names, **node_tag_names, **edge_tag_names, **face_tag_names, **cell_tag_names;
  real_t** field_data;
  polytope_tessellation_t* tess = polytope_tessellation_new(3);
  polytope_read_silo_with_tags(tess, &num_fields, &field_names, &field_data,
                               &num_node_tags, &node_tag_names, &node_tag_sizes, &node_tags,
                               &num_edge_tags, &edge_tag_names, &edge_tag_sizes, &edge_tags,
                               &num_face_tags, &face_tag_names, &face_tag_sizes, &face_tags,
                               &num_cell_tags, &cell_tag_names, &cell_tag_sizes, &cell_tags,
                               file_prefix, directory, cycle, time,
                               comm, num_files, mpi_tag);

  // Convert the data.

  // Mesh.
  *mesh = mesh_from_tessellation(tess);
  polytope_tessellation_free(tess);

  // Fields.
  *fields = string_ptr_unordered_map_new();
  for (int i = 0; i < num_fields; ++i)
    string_ptr_unordered_map_insert_with_kv_dtor(*fields, field_names[i], field_data[i], field_map_kv_dtor);
  free(field_data);
  free(field_names);

  // Tags.
  for (int i = 0; i < num_node_tags; ++i)
  {
    int* tag = mesh_create_tag((*mesh)->node_tags, node_tag_names[i], node_tag_sizes[i]);
    memcpy(tag, node_tags[i], sizeof(int) * node_tag_sizes[i]);
    free(node_tag_names[i]);
    free(node_tags[i]);
  }
  free(node_tag_names);
  free(node_tag_sizes);
  free(node_tags);
  for (int i = 0; i < num_edge_tags; ++i)
  {
    int* tag = mesh_create_tag((*mesh)->edge_tags, edge_tag_names[i], edge_tag_sizes[i]);
    memcpy(tag, edge_tags[i], sizeof(int) * edge_tag_sizes[i]);
    free(edge_tag_names[i]);
    free(edge_tags[i]);
  }
  free(edge_tag_names);
  free(edge_tag_sizes);
  free(edge_tags);
  for (int i = 0; i < num_face_tags; ++i)
  {
    int* tag = mesh_create_tag((*mesh)->face_tags, face_tag_names[i], face_tag_sizes[i]);
    memcpy(tag, face_tags[i], sizeof(int) * face_tag_sizes[i]);
    free(face_tag_names[i]);
    free(face_tags[i]);
  }
  free(face_tag_names);
  free(face_tag_sizes);
  free(face_tags);
  for (int i = 0; i < num_cell_tags; ++i)
  {
    int* tag = mesh_create_tag((*mesh)->cell_tags, cell_tag_names[i], cell_tag_sizes[i]);
    memcpy(tag, cell_tags[i], sizeof(int) * cell_tag_sizes[i]);
    free(cell_tag_names[i]);
    free(cell_tags[i]);
  }
  free(cell_tag_names);
  free(cell_tag_sizes);
  free(cell_tags);
}

