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
                    double* time,
                    MPI_Comm comm,
                    int num_files,
                    int mpi_tag)
{
  // Read the stuff from the file.
  int num_fields, num_node_tags, num_edge_tags, num_face_tags, num_cell_tags; 
  int *node_tag_sizes, *edge_tag_sizes, *face_tag_sizes, *cell_tag_sizes; 
  int **node_tags, **edge_tags, **face_tags, **cell_tags;
  char **field_names, **node_tag_names, **edge_tag_names, **face_tag_names, **cell_tag_names;
  double** field_data;
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

