// Copyright (c) 2012-2015, Jeffrey N. Johnson
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

#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"
#include "core/write_silo.h"
#include "core/array.h"

#if POLYMEC_HAVE_MPI
#include "mpi.h"
#include "pmpio.h"

static void* pmpio_create_file(const char* filename,
                               const char* dir_name,
                               void* user_data)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBMkDir(file, dir_name);
  DBSetDir(file, dir_name);
  return (void*)file;
}

static void* pmpio_open_file(const char* filename, 
                             const char* dir_name,
                             PMPIO_iomode_t iomode, 
                             void* user_data)
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

static void pmpio_close_file(void* file, void* user_data)
{
  DBClose((DBfile*)file);
}
#endif

static void write_tags_to_file(tagger_t* tagger, const char* tag_list_name, DBfile* file)
{
  // Pack the tags into a compound array.
  int_array_t* elem_lengths = int_array_new();
  string_array_t* elem_names = string_array_new();
  int_array_t* tag_data = int_array_new();

  int pos = 0, *tag, tag_size;
  char* tag_name;
  while (mesh_next_tag(tagger, &pos, &tag_name, &tag, &tag_size))
  {
    int_array_append(elem_lengths, tag_size);
    string_array_append(elem_names, tag_name);
    for (int i = 0; i < tag_size; ++i)
      int_array_append(tag_data, tag[i]);
  }

  // Write the compound array.
  if (elem_names->size > 0)
  {
    DBPutCompoundarray(file, tag_list_name, elem_names->data, elem_lengths->data,
                       elem_names->size, tag_data->data, tag_data->size, DB_INT, 0);
  }

  // Clean up.
  int_array_free(elem_lengths);
  string_array_free(elem_names);
  int_array_free(tag_data);
}

void write_silo_mesh(MPI_Comm comm,
                     const char* file_prefix,
                     const char* directory,
                     int cycle,
                     int num_files,
                     int mpi_tag,
                     mesh_t* mesh,
                     string_ptr_unordered_map_t* fields,
                     real_t time)
{
  // Strip .silo off of the prefix if it's there.
  char prefix[FILENAME_MAX];
  strncpy(prefix, file_prefix, FILENAME_MAX);
  char* suffix = strstr(prefix, ".silo");
  if (suffix != NULL)
    suffix[0] = '\0';

  // Open a file in Silo/HDF5 format for writing.
  char filename[FILENAME_MAX];
#if POLYMEC_HAVE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (num_files == -1)
    num_files = nproc;
  ASSERT(num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  char master_dir_name[FILENAME_MAX];
  if (strlen(directory) == 0)
    snprintf(master_dir_name, FILENAME_MAX, "%s-%d", prefix, nproc);
  else
    strncpy(master_dir_name, directory, FILENAME_MAX);
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
                                    pmpio_create_file, pmpio_open_file, 
                                    pmpio_close_file, 0);
  int group_rank = PMPIO_GroupRank(baton, rank);
  int rank_in_group = PMPIO_RankInGroup(baton, rank);

  // Create a subdirectory for each group.
  char group_dir_name[FILENAME_MAX];
  snprintf(group_dir_name, FILENAME_MAX, "%s/%d", master_dir_name, group_rank);
  if (rank_in_group == 0)
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
    snprintf(filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, prefix, cycle);
  else
    snprintf(filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, prefix);

  char dir_name[FILENAME_MAX];
  snprintf(dir_name, FILENAME_MAX, "domain_%d", rank_in_group);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dir_name);
#else
  char dir_name[FILENAME_MAX];
  if (strlen(directory) == 0)
    strcpy(dir_name, ".");
  else
    strcpy(dir_name, directory);

  if (cycle >= 0)
    snprintf(filename, FILENAME_MAX, "%s/%s-%d.silo", dir_name, prefix, cycle);
  else
    snprintf(filename, FILENAME_MAX, "%s/%s.silo", dir_name, prefix);

  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file, "/");
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  double dtime = (double)time;
  if (cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (dtime != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
  double* x = polymec_malloc(sizeof(double) * num_nodes);
  double* y = polymec_malloc(sizeof(double) * num_nodes);
  double* z = polymec_malloc(sizeof(double) * num_nodes);
  for (int i = 0; i < num_nodes; ++i)
  {
    x[i] = (double)mesh->nodes[i].x;
    y[i] = (double)mesh->nodes[i].y;
    z[i] = (double)mesh->nodes[i].z;
  }
  double* coords[3];
  coords[0] = x;
  coords[1] = y;
  coords[2] = z;

  // The polyhedral zone list is referred to in the options list.
  DBAddOption(optlist, DBOPT_PHZONELIST, (char*)"mesh_zonelist");

  // Write out the 3D polyhedral mesh.
  int num_cells = mesh->num_cells;
  int num_ghost_cells = mesh->num_ghost_cells;
  DBPutUcdmesh(file, (char*)"mesh", 3, coordnames, coords,
               num_nodes, num_cells + num_ghost_cells, 0, 0,
               DB_DOUBLE, optlist);

  // Partial cleanup.
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Construct the silo face-node info.  We rely on the mesh having
  // the faces nodes arranged counter-clockwise around the face.
  int num_faces = mesh->num_faces;
  int* face_node_counts = polymec_malloc(sizeof(int) * num_faces);
  char* ext_faces = polymec_malloc(sizeof(char) * num_faces);
  for (int i = 0; i < num_faces; ++i)
  {
    face_node_counts[i] = mesh->face_node_offsets[i+1] - mesh->face_node_offsets[i];
    if (mesh->face_cells[2*i+1] == -1)
      ext_faces[i] = 0x1;
    else
      ext_faces[i] = 0x0;
  }

  // Construct the silo cell-face info.  Silo uses the same 1's complement
  // convention we use for indicating face orientation, so we can
  // simply copy our faces.
  int* cell_face_counts = polymec_malloc(sizeof(int) * (num_cells + num_ghost_cells));
  memset(cell_face_counts, 0, sizeof(int) * (num_cells + num_ghost_cells));
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  DBPutPHZonelist(file, (char*)"mesh_zonelist", 
                  num_faces, face_node_counts,
                  mesh->face_node_offsets[num_faces], mesh->face_nodes,
                  ext_faces, num_cells + num_ghost_cells, cell_face_counts,
                  mesh->cell_face_offsets[num_cells], mesh->cell_faces,
                  0, 0, num_cells-1, optlist);

  // Partial cleanup.
  polymec_free(face_node_counts);
  polymec_free(ext_faces);
  polymec_free(cell_face_counts);

#if 0
  // Write out the cell-face connectivity data.
  int* conn = polymec_malloc(sizeof(int) * num_cells);
  int elem_lengths[3];
  char* elem_names[3];
  for (int c = 0; c < num_cells; ++c)
    conn[c] = mesh.cells[c].size();
  for (int c = 0; c < num_cells; ++c)
  {
    for (int f = 0; f < mesh.cells[c].size(); ++f) {
      int j = mesh.cells[c][f];
      conn.push_back(j < 0 ? ~j : j);
    }
  }
  for (int f = 0; f < mesh.faceCells.size(); ++f)
  {
    conn.push_back(mesh.faceCells[f][0]);
    conn.push_back(mesh.faceCells[f][0]);
  }
  elem_names[0] = strDup("ncellfaces");
  elem_lengths[0] = num_cells;
  elem_names[2] = strDup("facecells");
  elem_lengths[2] = conn.size() - 2*mesh.faces.size();
  elem_names[1] = strDup("cellfaces");
  elem_lengths[1] = conn.size() - elem_lengths[2] - elem_lengths[0];
  DBPutCompoundarray(file, "conn", elem_names, elem_lengths, 3, 
                     (void*)&conn[0], conn.size(), DB_INT, 0);
  polymec_free(elem_names[0]);
  polymec_free(elem_names[1]);
  polymec_free(elem_names[2]);
#endif

  // Write out tag information.
  write_tags_to_file(mesh->node_tags, "node_tags", file);
  write_tags_to_file(mesh->edge_tags, "edge_tags", file);
  write_tags_to_file(mesh->face_tags, "face_tags", file);
  write_tags_to_file(mesh->cell_tags, "cell_tags", file);

  // Write out the cell-centered field data.
  if (fields != NULL)
  {
    int pos = 0;
    char* field_name;
    void* item;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
    {
      DBPutUcdvar1(file, 
                   field_name,
                   (char*)"mesh",
                   item,
                   mesh->num_cells,
                   0, 
                   0,
                   DB_DOUBLE,
                   DB_ZONECENT,
                   optlist);
    }
  }

  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // Write the multi-block objects to the file if needed.
  int num_chunks = nproc / num_files;
  if (rank_in_group == 0)
  {
    char* mesh_names[num_chunks];
    int mesh_types[num_chunks];
    int var_types[num_chunks];
    char* var_names[num_chunks][fields->size];
    for (int i = 0; i < num_chunks; ++i)
    {
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "domain_%d/mesh", i);
      mesh_names[i] = string_dup(mesh_name);
      mesh_types[i] = DB_UCDMESH;
      var_types[i] = DB_UCDVAR;

      // Field data.
      int pos = 0, j = 0;
      char* field_name;
      void* item;
      while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
      {
        char var_name[1024];
        snprintf(var_name, 1024, "domain_%d/%s", i, field_name);
        var_names[j++][i] = string_dup(var_name);
      }
    }

    // Stick cycle and time in there if needed.
    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = (double)time;
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "mesh", num_chunks, &mesh_names[0], 
                   &mesh_types[0], optlist);
    {
      int pos = 0, j = 0;
      char* field_name;
      void* item;
      while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
        DBPutMultivar(file, field_name, num_chunks, var_names[j++], var_types, optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      polymec_free(mesh_names[i]);
    for (int f = 0; f < fields->size; ++f)
      for (int i = 0; i < num_chunks; ++i)
        polymec_free(var_names[f][i]);
  }

  // Write the file.
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_Finish(baton);

  // Finally, write the uber-master file.
  if (rank == 0)
  {
    char master_file_name[FILENAME_MAX];
    if (cycle >= 0)
      snprintf(master_file_name, FILENAME_MAX, "%s-%d/%s-%d.silo", prefix, nproc, prefix, cycle);
    else
      snprintf(master_file_name, FILENAME_MAX, "%s-%d/%s.silo", prefix, nproc, prefix);
    int driver = DB_HDF5;
    DBfile* file = DBCreate(master_file_name, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    char* mesh_names[num_files*num_chunks];
    int mesh_types[num_files*num_chunks];
    char* var_names[num_files*num_chunks][fields->size];
    int var_types[num_files*num_chunks];
    for (int i = 0; i < num_files; ++i)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        // Mesh.
        char mesh_name[1024];
        if (cycle >= 0)
          snprintf(mesh_name, 1024, "%d/%s-%d.silo:/domain_%d/mesh", i, prefix, cycle, c);
        else
          snprintf(mesh_name, 1024, "%d/%s.silo:/domain_%d/mesh", i, prefix, c);
        mesh_names[i*num_chunks+c] = string_dup(mesh_name);
        var_types[i*num_chunks+c] = DB_UCDVAR;

        // Field data.
        int pos = 0, j = 0;
        char* field_name;
        void* item;
        while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
        {
          char var_name[1024];
          if (cycle >= 0)
            snprintf(var_name, 1024, "%d/%s-%d.silo:/domain_%d/%s", i, prefix, cycle, c, field_name);
          else
            snprintf(var_name, 1024, "%d/%s.silo:/domain_%d/%s", i, prefix, c, field_name);
          var_names[j++][i*num_chunks+c] = string_dup(var_name);
        }
      }
    }

    DBoptlist* optlist = DBMakeOptlist(10);
    double dtime = (double)time;
    if (cycle >= 0)
      DBAddOption(optlist, DBOPT_CYCLE, &cycle);
    if (dtime != -FLT_MAX)
      DBAddOption(optlist, DBOPT_DTIME, &dtime);

    // Write the multimesh and variable data, and close the file.
    DBPutMultimesh(file, "mesh", num_files*num_chunks, &mesh_names[0], 
                   &mesh_types[0], optlist);
    {
      int pos = 0, j = 0;
      char* field_name;
      void* item;
      while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
        DBPutMultivar(file, field_name, num_files*num_chunks, var_names[j++], var_types, optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_files*num_chunks; ++i)
      polymec_free(mesh_names[i]);
    for (int f = 0; f < fields->size; ++f)
      for (int i = 0; i < num_files*num_chunks; ++i)
        polymec_free(var_names[f][i]);
  }
#else
  // Write the file.
  DBClose(file);
#endif
}

void write_silo_points(MPI_Comm comm,
                       const char* file_prefix,
                       const char* directory,
                       int cycle,
                       int num_files,
                       int mpi_tag,
                       point_t* points,
                       int num_points,
                       string_ptr_unordered_map_t* fields,
                       real_t time)
{
  // Floating point data type.
  int data_type = (sizeof(real_t) == sizeof(double)) ? DB_DOUBLE : DB_FLOAT;

  // Strip .silo off of the prefix if it's there.
  char prefix[FILENAME_MAX];
  strncpy(prefix, file_prefix, FILENAME_MAX);
  char* suffix = strstr(prefix, ".silo");
  if (suffix != NULL)
    suffix[0] = '\0';

  // Open a file in Silo/HDF5 format for writing.
  char filename[FILENAME_MAX];
#if POLYMEC_HAVE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (num_files == -1)
    num_files = nproc;
  ASSERT(num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  char master_dir_name[FILENAME_MAX];
  if (strlen(directory) == 0)
  {
    char dir_name[FILENAME_MAX];
    snprintf(dir_name, FILENAME_MAX, "%s-%d", prefix, nproc);
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
  int group_rank = PMPIO_GroupRank(baton, rank);
  int rank_in_group = PMPIO_RankInGroup(baton, rank);

  // Create a subdirectory for each group.
  char group_dir_name[FILENAME_MAX];
  snprintf(group_dir_name, FILENAME_MAX, "%s/%d", master_dir_name, group_rank);
  if (rank_in_group == 0)
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
    snprintf(filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, prefix, cycle);
  else
    snprintf(filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, prefix);

  char dir_name[FILENAME_MAX];
  snprintf(dir_name, FILENAME_MAX, "domain_%d", rank_in_group);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dir_name);
#else
  char dir_name[FILENAME_MAX];
  if (strlen(directory) == 0)
    strcpy(dir_name, ".");
  else
    strcpy(dir_name, directory);

  if (cycle >= 0)
    snprintf(filename, FILENAME_MAX, "%s/%s-%d.silo", dir_name, prefix, cycle);
  else
    snprintf(filename, FILENAME_MAX, "%s/%s.silo", dir_name, prefix);

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

  // Point coordinates.
  real_t* x = polymec_malloc(sizeof(real_t) * num_points);
  real_t* y = polymec_malloc(sizeof(real_t) * num_points);
  real_t* z = polymec_malloc(sizeof(real_t) * num_points);
  for (int i = 0; i < num_points; ++i)
  {
    x[i] = points[i].x;
    y[i] = points[i].y;
    z[i] = points[i].z;
  }
  real_t* coords[3];
  coords[0] = &(x[0]);
  coords[1] = &(y[0]);
  coords[2] = &(z[0]);

  // Write out the point mesh.
  DBPutPointmesh(file, (char*)"points", 3, coords, num_points, data_type, optlist); 
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Write out the point field data.
  if (fields != NULL)
  {
    int pos = 0;
    char* field_name; 
    real_t* field_data;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, (void**)&field_data))
    {
      real_t* vars[1] = {field_data}; 
      DBPutPointvar(file, field_name, "points", 1, vars, num_points, data_type, optlist);
    }
  }

  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // Write the multi-block objects to the file if needed.
  int num_chunks = nproc / num_files;
  int num_fields = fields->size;
  if (rank_in_group == 0)
  {
    char** mesh_names = polymec_malloc(sizeof(char*) * num_chunks);
    int* mesh_types = polymec_malloc(sizeof(int) * num_chunks);
    char*** var_names = polymec_malloc(sizeof(char**) * num_fields);
    for (int f = 0; f < num_fields; ++f)
      var_names[f] = polymec_malloc(sizeof(char*) * num_chunks);
    int* var_types = polymec_malloc(sizeof(int) * num_chunks);
    for (int i = 0; i < num_chunks; ++i)
    {
      mesh_types[i] = DB_POINTMESH;
      var_types[i] = DB_UCDVAR;
      // Mesh.
      char mesh_name[1024];
      snprintf(mesh_name, 1024, "domain_%d/points", i);
      mesh_names[i] = strdup(mesh_name);
      int pos = 0, field_index = 0;
      char* field_name;
      void* field_data;
      while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field_data))
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
    {
      double dtime = time;
      DBAddOption(optlist, DBOPT_DTIME, &dtime);
    }

    // Write the point mesh and variable data.
    DBSetDir(file, "/");
    DBPutMultimesh(file, "points", num_chunks, mesh_names, 
                   mesh_types, optlist);
    int pos = 0, field_index = 0;
    char* field_name;
    void* field_data;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field_data))
    {
      DBPutMultivar(file, field_name, num_chunks, var_names[field_index++], 
                    var_types, optlist);
    }

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      polymec_free(mesh_names[i]);
    polymec_free(mesh_names);
    polymec_free(mesh_types);
    for (int f = 0; f < num_fields; ++f)
    {
      for (int i = 0; i < num_chunks; ++i)
        polymec_free(var_names[f][i]);
      polymec_free(var_names[f]);
    }
    polymec_free(var_names);
    polymec_free(var_types);
  }

  // Write the file.
  PMPIO_HandOffBaton(baton, (void*)file);
  PMPIO_Finish(baton);

  // Finally, write the uber-master file.
  if (rank == 0)
  {
    char master_file_name[FILENAME_MAX];
    if (cycle >= 0)
      snprintf(master_file_name, FILENAME_MAX, "%s-%d/%s-%d.silo", prefix, nproc, prefix, cycle);
    else
      snprintf(master_file_name, FILENAME_MAX, "%s-%d/%s.silo", prefix, nproc, prefix);
    int driver = DB_HDF5;
    DBfile* file = DBCreate(master_file_name, DB_CLOBBER, DB_LOCAL, "Master file", driver);

    char** mesh_names = polymec_malloc(sizeof(char*) * num_files*num_chunks);
    int* mesh_types = polymec_malloc(sizeof(int) * num_files*num_chunks);
    char*** var_names = polymec_malloc(sizeof(char**) * num_fields);
    for (int f = 0; f < num_fields; ++f)
      var_names[f] = polymec_malloc(sizeof(char*) * num_files*num_chunks);
    int* var_types = polymec_malloc(sizeof(int) * num_files*num_chunks);
    for (int i = 0; i < num_files; ++i)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        mesh_types[num_chunks*i+c] = DB_POINTMESH;
        var_types[num_chunks*i+c] = DB_UCDVAR;

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
        void* field_data;
        while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field_data))
        {
          char var_name[1024];
          snprintf(var_name, 1024, "%d/domain_%d/%s", i, c, field_name);
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
    void* field_data;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field_data))
    {
      DBPutMultivar(file, field_name, num_chunks, var_names[field_index++], 
                    var_types, optlist);
    }
    DBClose(file);

    // Clean up.
    DBFreeOptlist(optlist);
    for (int i = 0; i < num_chunks; ++i)
      polymec_free(mesh_names[i]);
    polymec_free(mesh_names);
    polymec_free(mesh_types);
    for (int f = 0; f < num_fields; ++f)
    {
      for (int i = 0; i < num_files*num_chunks; ++i)
        polymec_free(var_names[f][i]);
      polymec_free(var_names[f]);
    }
    polymec_free(var_names);
    polymec_free(var_types);
  }
#else
  // Write the file.
  DBClose(file);
#endif
}


