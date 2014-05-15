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
#include "core/silo_file.h"
#include "core/array.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define SILO_FLOAT_TYPE DB_DOUBLE
#else
#define SILO_FLOAT_TYPE DB_FLOAT
#endif

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

// Object representing data in a "multi-mesh" or "multi-variable."
typedef struct
{
  char* mesh_name;
  int mesh_type;
  int var_type;
  char** field_names;
  int num_fields;
} multiobject_t;

// Constructors for various multi-objects.
static multiobject_t* multiobject_new(const char* mesh_name,
                                      int mesh_type,
                                      int var_type,
                                      string_ptr_unordered_map_t* fields)
{
  multiobject_t* obj = malloc(sizeof(multiobject_t));
  obj->mesh_name = string_dup(mesh_name);
  obj->mesh_type = mesh_type;
  obj->var_type = var_type;
  obj->num_fields = fields->size;
  obj->fields = malloc(sizeof(char*) * obj->num_fields);
  int pos = 0, i = 0;
  char* field_name;
  void* field;
  while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field))
    obj->fields[i++] = string_dup(field_name);
  return obj;
}

#endif

struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX], dir_name[FILENAME_MAX], filename[FILENAME_MAX];
  int cycle;
  real_t time;

#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  char master_dir_name[FILENAME_MAX], group_dir_name[FILENAME_MAX];
  PMPIO_baton_t* baton;
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;
  string_ptr_unordered_map_t* multiobjects;
#endif
};

#if POLYMEC_HAVE_MPI
static void write_multivars_to_file(silo_file_t* file)
{
  if (file->rank_in_group != 0) return;

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (file->time != -FLT_MAX)
  {
    double t = file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  int num_chunks = file->nproc / file->num_files;
  int pos = 0;
  multiobject_t* obj;
  while (string_ptr_unordered_map_next(file->multiobjects, &pos, &obj))
  {

    // Mesh and fields.
    char mesh_names[num_chunks][FILENAME_MAX];
    int mesh_types[num_chunks];
    char var_names[obj->num_fields][num_chunks][FILENAME_MAX];
    int var_types[num_chunks];
    for (int i = 0; i < num_chunks; ++i)
    {
      mesh_types[i] = obj->mesh_type;
      var_types[i] = obj->var_type;

      // Mesh.
      snprintf(mesh_names[i], FILENAME_MAX, "domain_%d/%s", i, obj->mesh_name);

      // Field names.
      for (int j = 0; j < obj->num_fields; ++j)
        snprintf(var_names[j][i], FILENAME_MAX, "domain_%d/%s_%s", i, obj->mesh_name, obj->field_names[j]);

      // Write the point mesh and variable data.
      DBSetDir(file, "/");
      DBPutMultimesh(file, "%s", obj­>mesh_name, num_chunks, mesh_names, 
          mesh_types, optlist);
      int field_index = 0;
      for (int j = 0; j < obj->num_fields; ++j, ++field_index)
      {
        DBPutMultivar(file, obj->field_names[j], num_chunks, var_names[field_index++], 
            var_types, optlist);
      }

    }
  }

  // Clean up.
  DBFreeOptlist(optlist);
}

static void write_master_file(silo_file_t* file)
{
  if (file->rank != 0) return;

  char master_file_name[FILENAME_MAX];
  snprintf(master_file_name, FILENAME_MAX, "%s-%d/%s.silo", file->prefix, file->nproc, file->prefix);
  int driver = DB_HDF5;
  DBfile* master = DBCreate(master_file_name, DB_CLOBBER, DB_LOCAL, "Master file", driver);

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &cycle);
  if (file->time != -FLT_MAX)
  {
    double t = file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  int num_files = file->num_files;
  int num_chunks = file->nproc / num_files;
  int pos = 0;
  multiobject_t* obj;
  while (string_ptr_unordered_map_next(file->multiobjects, &pos, &obj))
  {
    // Mesh and fields.
    char mesh_names[file->num_files*num_chunks][FILENAME_MAX];
    int mesh_types[file->num_files*num_chunks];
    char var_names[obj->num_fields][file->num_files*num_chunks][FILENAME_MAX];
    int var_types[num_files*num_chunks];
    for (int j = 0; j < file->num_files; ++j)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        mesh_types[num_chunks*j+c] = obj->mesh_type;
        var_types[num_chunks*j+c] = obj->var_type;

        // Mesh.
        snprintf(mesh_names[num_chunks*j+c], FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, prefix, c, obj->mesh_name);

        // Field names.
        for (int k = 0; k < obj->num_fields; ++k)
          snprintf(var_names[k][num_chunks*j+c], FILENAME_MAX, "%d/%s.silo:/domain_%d/%s_%s", j, file->prefix, c, file->meshes->data[i], fields->data[k]);
      }

      // Write the multimesh and variable data.
      DBPutMultimesh(file, "mesh", file->num_files*num_chunks, &mesh_names[0], 
          &mesh_types[0], optlist);
      int pos = 0, field_index = 0;
      char* field_name;
      void* field_data;
      while (string_ptr_unordered_map_next(fields, &pos, &field_name, &field_data))
      {
        DBPutMultivar(file, field_name, num_chunks, var_names[field_index++], 
            var_types, optlist);
      }
    }
  }

  DBFreeOptlist(optlist);
  DBClose(file);
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

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int num_files,
                            int mpi_tag)
{
  silo_file_t* file = malloc(sizeof(silo_file_t));
  file->cycle = -1;
  file->time = -FLT_MAX;

  // Strip .silo off of the prefix if it's there.
  {
    char prefix[FILENAME_MAX];
    strncpy(prefix, file_prefix, FILENAME_MAX);
    char* suffix = strstr(prefix, ".silo");
    if (suffix != NULL)
      suffix[0] = '\0';
    strcpy(file->prefix, prefix);
  }

#if POLYMEC_HAVE_MPI
  MPI_Comm_size(file->comm, &file->nproc);
  MPI_Comm_rank(file->comm, &file->rank);
  if (file->num_files == -1)
    file->num_files = nproc;
  else
    file->num_files = num_files;
  ASSERT(file->num_files <= nproc);

  // We put the entire data set into a directory named after the 
  // prefix, and every process gets its own subdirectory therein.

  // Create the master directory if we need to.
  if (strlen(directory) == 0)
    snprintf(file->master_dir_name, FILENAME_MAX, "%s-%d", file->prefix, file->nproc);
  else
    strncpy(file->master_dir_name, directory, FILENAME_MAX);
  if (rank == 0)
  {
    DIR* master_dir = opendir(file->master_dir_name);
    if (master_dir == NULL)
      mkdir(file->master_dir_name, S_IRWXU | S_IRWXG);
    else
      closedir(file->master_dir);
    MPI_Barrier(file->comm);
  }
  else
    MPI_Barrier(file->comm);

  // Initialize poor man's I/O and figure out group ranks.
  file->baton = PMPIO_Init(file->num_files, PMPIO_WRITE, file->comm, file->mpi_tag, 
                           pmpio_create_file, pmpio_open_file, 
                           pmpio_close_file, 0);
  file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
  file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

  // Create a subdirectory for each group.
  snprintf(file->group_dir_name, FILENAME_MAX, "%s/%d", file->master_dir_name, 
           file->group_rank);
  if (file->rank_in_group == 0)
  {
    DIR* group_dir = opendir(file->group_dir_name);
    if (group_dir == 0)
      mkdir(file->group_dir_name, S_IRWXU | S_IRWXG);
    else
      closedir(group_dir);
    MPI_Barrier(file->comm);
  }
  else
    MPI_Barrier(file->comm);

  // Determine a file name.
  snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->group_dir_name, prefix);

  snprintf(file->dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
  file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, file->dir_name);

  file->multiobjects = string_ptr_unordered_map_new();
#else
  if (strlen(directory) == 0)
    strncpy(file->dir_name, ".", FILENAME_MAX);
  else
    strncpy(file->dir_name, directory, FILENAME_MAX);

  snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->dir_name, file->prefix);

  int driver = DB_HDF5;
  file->dbfile = DBCreate(file->filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file->dbfile, "/");
#endif

  return file;
}

void silo_file_close(silo_file_t* file)
{
#if POLYMEC_HAVE_MPI
  // Write out the file itself.
  PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
  PMPIO_Finish(file->baton);

  // Write the uber-master file containing any multiobjects if need be.
  write_master_file(file);
#else
  // Write the file.
  DBClose(file->dbfile);
#endif
}

void silo_file_set_cycle(silo_file_t* file, int cycle)
{
  ASSERT(cycle >= 0);
  file->cycle = cycle;
}

void silo_file_set_time(silo_file_t* file, real_t time)
{
  ASSERT(time > -FLT_MAX);
  file->time = time;
}

void silo_file_add_mesh(silo_file_t* file,
                        const char* mesh_name,
                        mesh_t* mesh,
                        string_ptr_unordered_map_t* fields)
{
#if POLYMEC_HAVE_MPI
  file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, file->dir_name);
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  double dtime = (double)file->time;
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (dtime != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &dtime);

  // This is optional for now, but we'll give it anyway.
  char *coordnames[3];
  coordnames[0] = (char*)"xcoords";
  coordnames[1] = (char*)"ycoords";
  coordnames[2] = (char*)"zcoords";

  // Node coordinates.
  int num_nodes = mesh->num_nodes;
  double* x = malloc(sizeof(double) * num_nodes);
  double* y = malloc(sizeof(double) * num_nodes);
  double* z = malloc(sizeof(double) * num_nodes);
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
  char zonelist_name[FILENAME_MAX];
  snprintf(zonelist_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBAddOption(optlist, DBOPT_PHZONELIST, zonelist_name);

  // Write out the 3D polyhedral mesh.
  int num_cells = mesh->num_cells;
  int num_ghost_cells = mesh->num_ghost_cells;
  DBPutUcdmesh(file->dbfile, (char*)mesh_name, 3, coordnames, coords,
               num_nodes, num_cells + num_ghost_cells, 0, 0,
               SILO_FLOAT_TYPE, optlist);

  // Partial cleanup.
  free(x);
  free(y);
  free(z);

  // Construct the silo face-node info.  We rely on the mesh having
  // the faces nodes arranged counter-clockwise around the face.
  int num_faces = mesh->num_faces;
  int* face_node_counts = malloc(sizeof(int) * num_faces);
  char* ext_faces = malloc(sizeof(char) * num_faces);
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
  int* cell_face_counts = malloc(sizeof(int) * (num_cells + num_ghost_cells));
  memset(cell_face_counts, 0, sizeof(int) * (num_cells + num_ghost_cells));
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  DBPutPHZonelist(file->dbfile, zonelist_name, 
                  num_faces, face_node_counts,
                  mesh->face_node_offsets[num_faces], mesh->face_nodes,
                  ext_faces, num_cells + num_ghost_cells, cell_face_counts,
                  mesh->cell_face_offsets[num_cells], mesh->cell_faces,
                  0, 0, num_cells-1, optlist);

  // Partial cleanup.
  free(face_node_counts);
  free(ext_faces);
  free(cell_face_counts);

#if 0
  // Write out the cell-face connectivity data.
  int* conn = malloc(sizeof(int) * num_cells);
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
  free(elem_names[0]);
  free(elem_names[1]);
  free(elem_names[2]);
#endif

  // Write out tag information.
  write_tags_to_file(mesh->node_tags, "node_tags", file->dbfile);
  write_tags_to_file(mesh->edge_tags, "edge_tags", file->dbfile);
  write_tags_to_file(mesh->face_tags, "face_tags", file->dbfile);
  write_tags_to_file(mesh->cell_tags, "cell_tags", file->dbfile);

  // Write out the cell-centered field data.
  if (fields != NULL)
  {
    int pos = 0;
    char* field_name;
    void* item;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, &item))
    {
      DBPutUcdvar1(file->dbfile, 
                   field_name,
                   (char*)mesh_name,
                   item,
                   mesh->num_cells,
                   0, 
                   0,
                   SILO_FLOAT_TYPE,
                   DB_ZONECENT,
                   optlist);
    }
  }

  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // Write multi-block objects to the file if needed.
  write_multivars_to_file(file);
  PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
#endif
}

void silo_file_add_point_mesh(silo_file_t* file,
                              const char* point_mesh_name,
                              point_t* points,
                              int num_points,
                              string_ptr_unordered_map_t* fields)
{
#if POLYMEC_HAVE_MPI
  file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, file->dir_name);
#endif

  // Add cycle/time metadata if needed.
  DBoptlist* optlist = DBMakeOptlist(10);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  double t = file->time;
  if (t != -FLT_MAX)
    DBAddOption(optlist, DBOPT_DTIME, &t);

  // Point coordinates.
  real_t* x = malloc(sizeof(real_t) * num_points);
  real_t* y = malloc(sizeof(real_t) * num_points);
  real_t* z = malloc(sizeof(real_t) * num_points);
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
  DBPutPointmesh(file->dbfile, (char*)point_mesh_name, 3, coords, num_points, SILO_FLOAT_TYPE, optlist); 
  free(x);
  free(y);
  free(z);

  // Write out the point field data.
  if (fields != NULL)
  {
    int pos = 0;
    char* field_name; 
    real_t* field_data;
    while (string_ptr_unordered_map_next(fields, &pos, &field_name, (void**)&field_data))
    {
      real_t* vars[1] = {field_data}; 
      DBPutPointvar(file->dbfile, field_name, point_mesh_name, 1, vars, num_points, SILO_FLOAT_TYPE, optlist);
    }
  }

  // Clean up.
  DBFreeOptlist(optlist);

#if POLYMEC_HAVE_MPI
  // Write multi-block objects to the file if needed.
  write_multivars_to_file(file);
  PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
#endif
}

