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
#include "polytope_c.h"
#include "mpi.h"
#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

// Create a polytope tessellation from our mesh. This tessellation 
// borrows memory from the mesh, so it shouldn't be freed.
static polytope_tessellation_t* tessellation_from_mesh(mesh_t* mesh)
{
#ifndef NDEBUG
  // Check the validity of our mesh.
  for (int i = 0; i < mesh->num_cells; ++i)
  {
    ASSERT((mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i]) > 0);
  }
  for (int i = 0; i < mesh->num_faces; ++i)
  {
    ASSERT((mesh->face_node_offsets[i+1] - mesh->face_node_offsets[i]) > 0);
  }
#endif

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

void write_silo_mesh(mesh_t* mesh,
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

  // Translate the fields into polytope lingo.
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

  // Now fetch tags from the mesh and translate them too.

  int num_node_tags = 0;
  char* tag_name;
  int *tag_indices, tag_size, offset;
  pos = 0;
  while (mesh_next_tag(mesh->node_tags, &pos, &tag_name, &tag_indices, &tag_size))
    ++num_node_tags;
  char** node_tag_names = malloc(sizeof(char*) * MAX(1, num_node_tags));
  int* node_tag_sizes = malloc(sizeof(int) * MAX(1, num_node_tags));
  int** node_tag_indices = malloc(sizeof(int*) * MAX(1, num_node_tags));
  pos = offset = 0;
  while (mesh_next_tag(mesh->node_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    node_tag_names[offset] = tag_name; // borrowed!
    node_tag_sizes[offset] = tag_size;
    node_tag_indices[offset++] = tag_indices; // borrowed!
  }

  int num_edge_tags = 0;
  pos = 0;
  while (mesh_next_tag(mesh->edge_tags, &pos, &tag_name, &tag_indices, &tag_size))
    ++num_edge_tags;
  char** edge_tag_names = malloc(sizeof(char*) * MAX(1, num_edge_tags));
  int* edge_tag_sizes = malloc(sizeof(int) * MAX(1, num_edge_tags));
  int** edge_tag_indices = malloc(sizeof(int*) * MAX(1, num_edge_tags));
  pos = offset = 0;
  while (mesh_next_tag(mesh->edge_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    edge_tag_names[offset] = tag_name; // borrowed!
    edge_tag_sizes[offset] = tag_size;
    edge_tag_indices[offset++] = tag_indices; // borrowed!
  }

  int num_face_tags = 0;
  pos = 0;
  while (mesh_next_tag(mesh->face_tags, &pos, &tag_name, &tag_indices, &tag_size))
    ++num_face_tags;
  char** face_tag_names = malloc(sizeof(char*) * MAX(1, num_face_tags));
  int* face_tag_sizes = malloc(sizeof(int) * MAX(1, num_face_tags));
  int** face_tag_indices = malloc(sizeof(int*) * MAX(1, num_face_tags));
  pos = offset = 0;
  while (mesh_next_tag(mesh->face_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    face_tag_names[offset] = tag_name; // borrowed!
    face_tag_sizes[offset] = tag_size;
    face_tag_indices[offset++] = tag_indices; // borrowed!
  }

  int num_cell_tags = 0;
  pos = 0;
  while (mesh_next_tag(mesh->cell_tags, &pos, &tag_name, &tag_indices, &tag_size))
    ++num_cell_tags;
  char** cell_tag_names = malloc(sizeof(char*) * MAX(1, num_cell_tags));
  int* cell_tag_sizes = malloc(sizeof(int) * MAX(1, num_cell_tags));
  int** cell_tag_indices = malloc(sizeof(int*) * MAX(1, num_cell_tags));
  pos = offset = 0;
  while (mesh_next_tag(mesh->cell_tags, &pos, &tag_name, &tag_indices, &tag_size))
  {
    cell_tag_names[offset] = tag_name; // borrowed!
    cell_tag_sizes[offset] = tag_size;
    cell_tag_indices[offset++] = tag_indices; // borrowed!
  }

  polytope_write_silo_with_tags(tess, 
                                num_node_fields, node_field_names, node_field_data,
                                num_node_tags, node_tag_names, node_tag_sizes, node_tag_indices,
                                num_edge_fields, edge_field_names, edge_field_data,
                                num_edge_tags, edge_tag_names, edge_tag_sizes, edge_tag_indices,
                                num_face_fields, face_field_names, face_field_data,
                                num_face_tags, face_tag_names, face_tag_sizes, face_tag_indices,
                                num_cell_fields, cell_field_names, cell_field_data,
                                num_cell_tags, cell_tag_names, cell_tag_sizes, cell_tag_indices,
                                file_prefix, directory, cycle, time,
                                comm, num_files, mpi_tag);

  // Clean up.
  free(cell_tag_indices);
  free(cell_tag_sizes);
  free(cell_tag_names);
  free(face_tag_indices);
  free(face_tag_sizes);
  free(face_tag_names);
  free(edge_tag_indices);
  free(edge_tag_sizes);
  free(edge_tag_names);
  free(node_tag_indices);
  free(node_tag_sizes);
  free(node_tag_names);
  polytope_tessellation_free(tess);
}

#if HAVE_MPI
static void* pmpio_create_file(const char* filename,
                               const char* dir_name,
                               void* userData)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dir_name);
  DBSetDir(file, dir_name);
  return (void*)file;
}

void*
pmpio_open_file(const char* filename, 
                const char* dir_name,
                PMPIO_iomode_t iomode, 
                void* userData)
{
  int driver = DB_HDF5;
  DBfile* file;
  if (iomode == PMPIO_WRITE)
  { 
    file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
    DBMkDir(file, dir_name);
    DBSetDir(file, dir_name);
  }
  else
  {
    file = DBOpen(filename, driver, DB_READ);
    DBSetDir(file, dir_name);
  }
  return (void*)file;
}

void
pmpio_close_file(void* file,
                void* userData)
{
  DBClose((DBfile*)file);
}
#endif

void write_silo_points(point_t* points,
                       int num_points,
                       string_ptr_unordered_map_t* point_fields,
                       const char* file_prefix,
                       const char* directory,
                       int cycle,
                       double time,
                       MPI_Comm comm,
                       int num_files,
                       int mpi_tag)
{
  // Strip .silo off of the prefix if it's there.
  char* prefix = string_dup(file_prefix);
  char* suffix = strstr(prefix, ".silo");
  if (suffix != NULL)
    suffix[0] = '\0';

  // Open a file in Silo/HDF5 format for writing.
  char filename[1024];
#if HAVE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (num_files == -1)
    num_files = nproc;
  ASSERT(num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  char master_dir_name[1024];
  if (strlen(directory) == 0)
  {
    char dir_name[1024];
    snprintf(dir_name, 1024, "%s-%d", prefix, nproc);
    strcpy(master_dir_name, dir_name);
  }
  else
    strcpy(master_dir_name, directory);
  if (rank == 0)
  {
    DIR* master_dir = opendir(master_dir_name);
    if (master_dir == 0)
      mkdir(master_dir_name, S_IRWXU | S_IRWXG);
    else
      closedir(master_dir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(num_files, PMPIO_WRITE, comm, mpi_tag, 
                                    &pmpio_create_file, 
                                    &pmpio_open_file, 
                                    &pmpio_close_file,
                                    0);
  int group_rank = PMPIO_group_rank(baton, rank);
  int rankInGroup = PMPIO_RankInGroup(baton, rank);

  // Create a subdirectory for each group.
  char group_dir_name[1024];
  snprintf(group_dir_name, 1024, "%s/%d", masterdir_name, group_rank);
  if (rankInGroup == 0)
  {
    DIR* group_dir = opendir(group_dir_name);
    if (group_dir == 0)
      mkdir((char*)group_dir_name, S_IRWXU | S_IRWXG);
    else
      closedir(group_dir);
    MPI_Barrier(comm);
  }
  else
  {
    MPI_Barrier(comm);
  }

  // Determine a file name.
  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", group_dir_name, prefix, cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", group_dir_name, prefix);

  char dir_name[1024];
  snprintf(dir_name, 1024, "domain_%d", rankInGroup);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dir_name);
#else
  char dir_name[1024];
  if (strlen(directory) == 0)
    strcpy(dir_name, ".");
  else
    strcpy(dir_name, directory);

  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", dir_name, prefix, cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", dir_name, prefix);

  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file, "/");
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  if (cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (time != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &time);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Point coordinates.
  double* x = malloc(sizeof(double) * num_points);
  double* y = malloc(sizeof(double) * num_points);
  double* z = malloc(sizeof(double) * num_points);
  for (int i = 0; i < num_points; ++i)
  {
    x[i] = points[i].x;
    y[i] = points[i].y;
    z[i] = points[i].z;
  }
  double* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Write out the point mesh.
  DBPutPointmesh(file, (char*)"points", 3, coords, num_points, DB_DOUBLE, optlist); 
  free(x);
  free(y);
  free(z);

  // Write out the point field data.
  int pos = 0;
  char* field_name; 
  double* field_data;
  while (string_ptr_unordered_map_next(point_fields, &pos, &field_name, (void**)&field_data))
  {
    double* vars[1] = {field_data}; 
    DBPutPointvar(file, field_name, "points", 1, vars, num_points, DB_DOUBLE, optlist);
  }

  // Clean up.
  DBFreeOptlist(optlist);

#if HAVE_MPI
  // Write the multi-block objects to the file if needed.
  int num_chunks = nproc / num_files;
  int num_fields = point_fields->size;
  if (rankInGroup == 0)
  {
    char** mesh_names = malloc(sizeof(char*) * num_chunks);
    int* mesh_types = malloc(sizeof(int) * num_chunks);
    char*** var_names = malloc(sizeof(char**) * num_fields);
    for (int f = 0; f < num_fields; ++f)
      var_names[f] = malloc(sizeof(char*) * num_chunks);
    int* var_types = malloc(sizeof(int) * num_chunks);
    for (int i = 0; i < num_chunks; ++i)
    {
      mesh_types[i] = DB_POINTMESH;
      var_types[i] = DB_UCDV;
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "domain_%d/points", i);
      mesh_names[i] = strdup(mesh_name);
      int pos = 0, field_index = 0;
      char* field_name;
      double* field_data;
      while (string_ptr_unordered_map_next(point_fields, &pos, &field_name, &field_data))
      {
        char var_name[1024];
        snprintf(var_name, 1024, "domain_%d/%s", i, field_name);
        var_names[field_index++][i] = string_dup(var_name);
      }
    }

    // Stick cycle and time in there if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (time != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &time);

    // Write the point mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "points", num_chunks, mesh_names, 
                   mesh_types, optlist);
    int pos = 0, field_index = 0;
    char* field_name;
    double* field_data;
    while (string_ptr_unordered_map_next(point_fields, &pos, &field_name, &field_data))
    {
      DBPutMultivar(file, field_name, num_chunks, var_names[field_index++], 
                    var_types, optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      free(mesh_names[i]);
    free(mesh_names);
    free(mesh_types);
    for (int f = 0; f < var_names.size(); ++f)
    {
      for (int i = 0; i < num_chunks; ++i)
        free(var_names[f][i]);
      free(var_names[f]);
    }
    free(var_names);
    free(var_types);
  }

  // Write the file.
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_Finish(baton);

  // Finally, write the uber-master file.
  if (rank == 0)
  {
    char master_file_name[1024];
    if (cycle >= 0)
      snprintf(master_file_name, 1024, "%s-%d/%s-%d.silo", prefix, nproc, prefix, cycle);
    else
      snprintf(master_file_name, 1024, "%s-%d/%s.silo", prefix, nproc, prefix);
    int driver = DB_HDF5;
    DBfile* file = DBCreate(master_file_name, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    vector<char*> mesh_names(num_files*num_chunks);
    vector<int> mesh_types(num_files*num_chunks, DB_UCDMESH);
    vector<vector<char*> > var_names(nodeFields.size() +
                                    edgeFields.size() +
                                    faceFields.size() +
                                    cellFields.size());
    vector<int> varTypes(num_files*num_chunks, DB_UCDVAR);
    char** mesh_names = malloc(sizeof(char*) * num_files*num_chunks);
    int* mesh_types = malloc(sizeof(int) * num_files*num_chunks);
    char*** var_names = malloc(sizeof(char**) * num_fields);
    for (int f = 0; f < num_fields; ++f)
      var_names[f] = malloc(sizeof(char*) * num_files*num_chunks);
    int* var_types = malloc(sizeof(int) * num_files*num_chunks);
    for (int i = 0; i < num_files; ++i)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        mesh_types[num_chunks*i+c] = DB_POINTMESH;
        var_types[num_chunks*i+c] = DB_UCDV;

        // Mesh.
        char mesh_name[1024];
        if (cycle >= 0)
          snprintf(mesh_name, 1024, "%d/%s-%d.silo:/domain_%d/points", i, prefix, cycle, c);
        else
          snprintf(mesh_name, 1024, "%d/%s.silo:/domain_%d/points", i, prefix, c);
        mesh_names[num_chunks*i+c] = string_dup(mesh_name);

        // Field names.
        int pos = 0, field_index = 0;
        char* field_name;
        double* field_data;
        while (string_ptr_unordered_map_next(point_fields, &pos, &field_name, &field_data))
        {
          char var_name[1024];
          snprintf(var_name, 1024, "%d/domain_%d/%s", i, field_name);
          if (cycle >= 0)
            snprintf(var_name, 1024, "%d/%s-%d.silo:/domain_%d/%s", i, prefix, cycle, c, field_name);
          else
            snprintf(var_name, 1024, "%d/%s.silo:/domain_%d/%s", i, prefix, c, field_name);
          var_names[field_index++][num_chunks*i+c] = string_dup(var_name);
        }
      }
    }

    DBoptlist* optlist = DBMakeOptlist(10);
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (time != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &time);

    // Write the multimesh and variable data, and close the file.
    DBPutMultimesh(file, "mesh", num_files*num_chunks, &mesh_names[0], 
                   &mesh_types[0], optlist);
    int pos = 0, field_index = 0;
    char* field_name;
    double* field_data;
    while (string_ptr_unordered_map_next(point_fields, &pos, &field_name, &field_data))
    {
      DBPutMultivar(file, field_name, num_chunks, var_names[field_index++], 
                    var_types, optlist);
    }
    DBClose(file);

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      free(mesh_names[i]);
    free(mesh_names);
    free(mesh_types);
    for (int f = 0; f < num_fields; ++f)
    {
      for (int i = 0; i < num_files*num_chunks; ++i)
        free(var_names[f][i]);
      free(var_names[f]);
    }
    free(var_names);
    free(var_types);
  }
#else
  // Write the file.
  DBClose(file);
#endif
}


