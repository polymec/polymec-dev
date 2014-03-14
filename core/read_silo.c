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

#include <sys/stat.h>
#include <dirent.h>
#include "silo.h"

#include "core/read_silo.h"
#include "core/slist.h"

#if POLYMEC_HAVE_MPI

#include "mpi.h"
#include "pmpio.h"

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

static void* pmpio_open_file(const char* filename, 
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

static void pmpio_close_file(void* file, void* userData)
{
  DBClose((DBfile*)file);
}

#endif

static void field_map_kv_dtor(char* key, void* value)
{
  free(key);
  free(value);
}

void read_silo_mesh(MPI_Comm comm,
                    const char* file_prefix,
                    const char* directory,
                    int cycle,
                    int num_files,
                    int mpi_tag,
                    mesh_t** mesh,
                    string_ptr_unordered_map_t* fields,
                    real_t* time)
{
  // Strip .silo off of the prefix if it's there.
  char prefix[strlen(file_prefix)+1];
  char* tok = strstr(file_prefix, ".silo");
  if (tok != NULL)
    strncpy(prefix, file_prefix, (size_t)(tok - file_prefix));
  else
    strcpy(prefix, file_prefix);

  // Open a file in Silo/HDF5 format for reading.
  char filename[1024];
#if HAVE_MPI
  int nproc = 1, rank = 0;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  if (num_files == -1)
    num_files = nproc;
  POLY_ASSERT(num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Check for the existence of the master directory.
  char master_dir_name[strlen(directory)+1];
  strcpy(master_dir_name, directory);
  if (strlen(master_dir_name) == 0)
  {
    char dir_name[1024];
    snprintf(dir_name, 1024, "%s-%d", prefix, nproc);
    strcpy(master_dir_name, dir_name);
  }
  DIR* master_dir = opendir(master_dir_name);
  if (master_dir == 0)
    polymec_error("Could not find the directory %s", master_dir_name);

  // Initialize poor man's I/O and figure out group ranks.
  PMPIO_baton_t* baton = PMPIO_Init(num_files, PMPIO_READ, comm, mpi_tag, 
                                    pmpio_create_file, pmpio_open_file, 
                                    pmpio_close_file, 0);
  int group_rank = PMPIO_group_rank(baton, rank);
  int rank_in_group = PMPIO_rank_in_group(baton, rank);

  // Figure out the subdirectory for this group.
  char group_dir_name[1024];
  snprintf(group_dir_name, 1024, "%s/%d", master_dir_name, group_rank);
  DIR* group_dir = opendir(group_dir_name);
  if (group_dir == 0)
    polymec_error("Could not find the directory %s", group_dir_name);

  // Determine the file name.
  if (cycle >= 0)
    snprintf(filename, 1024, "%s/%s-%d.silo", group_dir_name, prefix, cycle);
  else
    snprintf(filename, 1024, "%s/%s.silo", group_dir_name, prefix);

  char dir_name[1024];
  snprintf(dir_name, 1024, "domain_%d", rank_in_group);
  DBfile* file = (DBfile*)PMPIO_WaitForBaton(baton, filename, dir_name);
  DBSetDir(file, dir_name);
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
  DBfile* file = DBOpen(filename, driver, DB_READ);
  DBSetDir(file, "/");
#endif

  // Retrieve the mesh. Note that we must deallocate the storage 
  // for this object after we're through!
  DBucdmesh* dbmesh = DBGetUcdmesh(file, "mesh");

  // Extract time.
  *time = dbmesh->dtime;

  // Initialize the mesh.
  int num_cells = dbmesh->zones->nzones;
  int num_ghost_cells = 0; // FIXME: dbmesh->zones->hi_offset - dbmesh->zones->lo_offset;
  int num_faces = dbmesh->faces->nfaces;
  int num_edges = 0; // To be constructed.
  int num_nodes = dbmesh->nnodes;
  mesh_t* mesh_obj = mesh_new(comm, num_cells, num_ghost_cells, num_faces, 
                              num_edges, num_nodes);

  // Node coordinates.
  for (int i = 0; i < dbmesh->nnodes; ++i)
  {
    mesh_obj->nodes[i].x = ((real_t*)(dbmesh->coords[i]))[0];
    mesh_obj->nodes[i].y = ((real_t*)(dbmesh->coords[i]))[1];
    mesh_obj->nodes[i].z = ((real_t*)(dbmesh->coords[i]))[2];
  }

  // Reconstruct the faces.
  mesh_obj->face_node_offsets[0] = 0;
  for (int f = 0; f < mesh_obj->num_faces; ++f)
    mesh_obj->face_node_offsets[f+1] = mesh_obj->face_node_offsets[f] + dbmesh->faces->shapesize[f];
  mesh_obj->face_nodes = ARENA_REALLOC(mesh_obj->arena, mesh_obj->face_nodes, sizeof(int)*mesh_obj->face_node_offsets[mesh_obj->num_faces], 0);
  memcpy(mesh_obj->face_nodes, dbmesh->faces->nodelist, sizeof(int)*mesh_obj->face_node_offsets[mesh_obj->num_faces]);

  // Reconstruct the cell-face connectivity.
  DBcompoundarray* conn = DBGetCompoundarray(file, "connectivity");
  if (conn == 0)
  {
    DBClose(file);
    polymec_error("Could not find cell-face connectivity in file %s.", filename);
  }
  // First element is the number of faces in each zone.
  // Second element is the list of face indices in each zone.
  // Third element is a pair of cells for each face.
  if ((conn->nelems != 3) || 
      (conn->elemlengths[0] != dbmesh->zones->nzones) || 
      (conn->elemlengths[2] != 2*dbmesh->faces->nfaces))
  {
    DBClose(file);
    polymec_error("Found invalid cell-face connectivity in file %s.", filename);
  }
  int* conn_data = conn->values;

  // Cell-face connectivity.
  mesh_obj->cell_face_offsets[0] = 0;
  for (int c = 0; c < dbmesh->zones->nzones; ++c)
    mesh_obj->cell_face_offsets[c+1] = mesh_obj->cell_face_offsets[c] + conn_data[c];
  int foffset = dbmesh->zones->nzones;
  mesh_obj->cell_faces = ARENA_REALLOC(mesh_obj->arena, mesh_obj->cell_faces, sizeof(int)*mesh_obj->cell_face_offsets[mesh_obj->num_cells], 0);
  memcpy(mesh_obj->cell_faces, &conn_data[foffset], sizeof(int) * mesh_obj->cell_face_offsets[mesh_obj->num_cells]);
  foffset += mesh_obj->cell_face_offsets[mesh_obj->num_cells];

  // Face-cell connectivity.
  memcpy(mesh_obj->face_cells, &conn_data[foffset], sizeof(int)*mesh_obj->num_faces);

  // Close the mesh portion of the file.
  DBFreeUcdmesh(dbmesh);
  DBFreeCompoundarray(conn);

  // Construct mesh edges from existing data.
  mesh_construct_edges(mesh_obj);
  
  // Read any tag data.
  DBcompoundarray* tags = DBGetCompoundarray(file, "node_tags");
  if (tags != NULL)
  {
    int offset = 0;
    for (int i = 0; i < tags->nelems; ++i)
    {
      int* tag = mesh_create_tag(mesh_obj->node_tags, tags->elemnames[i], tags->elemlengths[i]);
      memcpy(tag, &tags->values[offset], sizeof(int) * tags->elemlengths[i]);
      offset += tags->elemlengths[i];
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "edge_tags");
  if (tags != NULL)
  {
    int offset = 0;
    for (int i = 0; i < tags->nelems; ++i)
    {
      int* tag = mesh_create_tag(mesh_obj->edge_tags, tags->elemnames[i], tags->elemlengths[i]);
      memcpy(tag, &tags->values[offset], sizeof(int) * tags->elemlengths[i]);
      offset += tags->elemlengths[i];
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "face_tags");
  if (tags != NULL)
  {
    int offset = 0;
    for (int i = 0; i < tags->nelems; ++i)
    {
      int* tag = mesh_create_tag(mesh_obj->face_tags, tags->elemnames[i], tags->elemlengths[i]);
      memcpy(tag, &tags->values[offset], sizeof(int) * tags->elemlengths[i]);
      offset += tags->elemlengths[i];
    }
    DBFreeCompoundarray(tags);
  }
  tags = DBGetCompoundarray(file, "cell_tags");
  if (tags != NULL)
  {
    int offset = 0;
    for (int i = 0; i < tags->nelems; ++i)
    {
      int* tag = mesh_create_tag(mesh_obj->cell_tags, tags->elemnames[i], tags->elemlengths[i]);
      memcpy(tag, &tags->values[offset], sizeof(int) * tags->elemlengths[i]);
      offset += tags->elemlengths[i];
    }
    DBFreeCompoundarray(tags);
  }

  // Retrieve the fields.
  int pos = 0;
  char* field_name;
  void* item;
  string_slist_t* field_names = string_slist_new();
  ptr_slist_t* new_field_data = ptr_slist_new();
  while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
  {
    ASSERT(item == NULL); // Should be unassociated, except for the name.
    DBucdvar* dbvar = DBGetUcdvar(file, field_name);
    if (dbvar == 0)
    {
      DBClose(file);
      polymec_error("Could not find field %s in file %s.", field_name, filename);
    }
    real_t* new_data = malloc(sizeof(real_t) * dbvar->nels);
    memcpy(new_data, dbvar->vals, sizeof(real_t) * dbvar->nels);
    string_slist_append(field_names, field_name);
    ptr_slist_append(new_field_data, new_data);

    // Clean up.
    DBFreeUcdvar(dbvar);
  }

  // Close the file
  DBClose(file);

  // Move things into place.
  *mesh = mesh_obj;
  string_slist_node_t* name_iter;
  ptr_slist_node_t* data_iter;
  for (name_iter = field_names->front, data_iter = new_field_data->front;
       name_iter != NULL;
       name_iter = name_iter->next, data_iter = data_iter->next)
  {
    string_ptr_unordered_map_insert_with_kv_dtor(fields, string_dup(name_iter->value), 
                                                 data_iter->value, field_map_kv_dtor);
  }
}

