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
#include "core/slist.h"
#include "core/array.h"

#if POLYMEC_HAVE_DOUBLE_PRECISION
#define SILO_FLOAT_TYPE DB_DOUBLE
#else
#define SILO_FLOAT_TYPE DB_FLOAT
#endif

#if POLYMEC_HAVE_MPI
#include "mpi.h"
#include "pmpio.h"

void silo_file_query(const char* file_prefix,
                     const char* directory,
                     int* num_files,
                     int* num_mpi_processes,
                     int* cycles,
                     int max_cycles_size)
{
  // No blank strings allowed for queries.
  ASSERT(strlen(file_prefix) > 0);
  ASSERT(strlen(directory) > 0);
  ASSERT(max_cycles_size >= 0);

  // Inspect the given directory's contents.
  if (!directory_exists(directory))
    polymec_error("Can't query non-existent directory: %s", directory);
  string_slist_t* files_in_dir = files_within_directory(directory);

  // Try to find a master file or single data file.
  char data_file[FILENAME_MAX];
  data_file[0] = '\0';
  bool found_cycles = false;
  {
    string_slist_node_t* node = files_in_dir->front;
    char path[FILENAME_MAX];
    snprintf(path, FILENAME_MAX, "%s/%s", directory, file_prefix); 
    if (strstr(node->value, path) && strstr(node->value, "silo"))
    {
      strncpy(data_file, node->value, FILENAME_MAX);
      snprintf(path, FILENAME_MAX, "%s/%s-", directory, file_prefix); 
      if (strstr(node->value, path))
        found_cycles = true;
    }
  }

  // Open up the file and see whether it's a master file or a serial data file.
  bool is_master = false;
  int driver = DB_HDF5;
  DBfile* file = DBOpen(data_file, driver, DB_READ);

  // What's in there?
  DBtoc* toc = DBGetToc(file); 
  if ((toc->nucdmesh == 0) && (toc->nptmesh == 0) && (toc->nmultimesh > 0))
    is_master = true;
  if (is_master)
  {
    // How many files are referenced in the master file?
    int fmax = -1;
    for (int f = 0; f < toc->nmultimesh; ++f)
    {
      char* p1 = toc->multimesh_names[f];
      char* p2 = strstr(toc->multimesh_names[f], "/");
      ASSERT(p2 != NULL);
      ASSERT(p2 > p1);
      char file_num[p2-p1+1];
      strncpy(file_num, p1, p2-p1);
      ASSERT(string_is_number(file_num));
      fmax = MAX(fmax, atoi(file_num));
    }
    *num_files = fmax + 1;

    // Grab one of the meshes and use its name as a template for identifying
    // the number of MPI processes used to construct the data set.
    {
      char* p1 = strstr(toc->multimesh_names[0], "domain_");
      ASSERT(p1 != NULL);
      char* p2 = p1 + strlen("domain_");
      char multimesh_prefix[p2-p1+1];
      strncpy(multimesh_prefix, toc->multimesh_names[0], p2-p1);
      int num_domains = 0;
      for (int f = 0; f < toc->nmultimesh; ++f)
      {
        char* s = strstr(toc->multimesh_names[f], multimesh_prefix);
        if (s != NULL)
          ++num_domains;
      }
      *num_mpi_processes = num_domains;
    }
  }
  else
  {
    // A single data file can only be written for a serial run.
    *num_files = 1;
    *num_mpi_processes = 1;
  }

  DBClose(file);

  // Search for available cycles.
  if ((cycles != NULL) && (max_cycles_size > 0) && found_cycles)
  {
    int num_cycles = 0;
    string_slist_node_t* node = files_in_dir->front;
    while (node != NULL)
    {
      char path[FILENAME_MAX];
      snprintf(path, FILENAME_MAX, "%s/%s-", directory, file_prefix); 
      char* p1 = strstr(node->value, path);
      char* p2 = strstr(node->value, ".silo");
      if ((p1 != NULL) && (p2 != NULL))
      {
        char* c = p1 + strlen(path);
        char num[p2 - c+1];
        strncpy(num, c, p2 - c);
        if (string_is_number(num))
          cycles[num_cycles++] = atoi(num);
        if (num_cycles == max_cycles_size)
          break;
      }
      node = node->next;
    }
  }

  // Clean up.
  string_slist_free(files_in_dir);
}

static void* pmpio_create_file(const char* filename,
                               const char* dir_name,
                               void* user_data)
{
  int driver = DB_HDF5;
  DBfile* file = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
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
    FILE* f = fopen(filename, "r");
    if (f == NULL)
      file = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
    else
    {
      fclose(f);
      file = DBOpen(filename, driver, DB_APPEND);
    }
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

// Object representing data in a multi-mesh.
typedef struct
{
  char* name;
  int type;
} multimesh_t;

static multimesh_t* multimesh_new(const char* mesh_name, int mesh_type)
{
  multimesh_t* mesh = polymec_malloc(sizeof(multimesh_t));
  mesh->name = string_dup(mesh_name);
  mesh->type = mesh_type;
  return mesh;
}

static void multimesh_free(multimesh_t* mesh)
{
  polymec_free(mesh->name);
  polymec_free(mesh);
}

// Object representing data in a multi-mesh.
typedef struct
{
  char* mesh_name;
  char* name;
  int type;
} multivar_t;

// Constructors for various multi-objects.
static multivar_t* multivar_new(const char* mesh_name,
                                const char* var_name,
                                int var_type)
{
  multivar_t* var = polymec_malloc(sizeof(multivar_t));
  var->mesh_name = string_dup(mesh_name);
  var->name = string_dup(var_name);
  var->type = var_type;
  return var;
}

static void multivar_free(multivar_t* var)
{
  polymec_free(var->mesh_name);
  polymec_free(var->name);
  polymec_free(var);
}

#endif

struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX], directory[FILENAME_MAX], filename[FILENAME_MAX];
  int cycle;
  real_t time;
  int mode; // Open for reading (DB_READ) or writing (DB_CLOBBER)? 

#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  PMPIO_baton_t* baton;
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;
  ptr_array_t* multimeshes;
  ptr_array_t* multivars;
#endif
};

#if POLYMEC_HAVE_MPI
static void write_multivars_to_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);

  if (file->rank_in_group != 0) return;

  int num_chunks = file->nproc / file->num_files;

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (file->time != -FLT_MAX)
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  // Write out multi meshes.
  for (int i = 0; i < file->multimeshes->size; ++i)
  {
    multimesh_t* mesh = file->multimeshes->data[i];

    char* mesh_names[num_chunks];
    int mesh_types[num_chunks];
    for (int j = 0; j < num_chunks; ++j)
    {
      char mesh_name[FILENAME_MAX];
      snprintf(mesh_name, FILENAME_MAX, "domain_%d/%s", j, mesh->name);
      mesh_names[j] = string_dup(mesh_name);
      mesh_types[j] = mesh->type;
    }

    // Write the point mesh and variable data.
    DBSetDir(file->dbfile, "/");
    DBPutMultimesh(file->dbfile, mesh->name, num_chunks, &mesh_names[0], 
                   mesh_types, optlist);

    // Clean up.
    for (int j = 0; j < num_chunks; ++j)
      polymec_free(mesh_names[j]);
  }

  // Multi variables.
  for (int i = 0; i < file->multivars->size; ++i)
  {
    multivar_t* var = file->multivars->data[i];

    // Fields and associated meshes.
    char* var_names[num_chunks];
    int var_types[num_chunks];
    for (int j = 0; j < num_chunks; ++j)
    {
      // Field name.
      char var_name[FILENAME_MAX];
      snprintf(var_name, FILENAME_MAX, "domain_%d/%s", j, var->name);
      var_names[j] = string_dup(var_name);
      var_types[j] = var->type;
    }

    // Write the variable data.
    DBSetDir(file->dbfile, "/");
    DBPutMultivar(file->dbfile, var->name, num_chunks, var_names, var_types, optlist);

    // Clean up.
    for (int j = 0; j < num_chunks; ++j)
      polymec_free(var_names[j]);
  }
  DBFreeOptlist(optlist);
}

static void write_master_file(silo_file_t* file)
{
  ASSERT(file->mode == DB_CLOBBER);

  // FIXME: Should change this to use Silo's name schemes for multi-block 
  // FIXME: objects when we start to Get Real Parallel.

  if (file->rank != 0) return;

  char master_file_name[FILENAME_MAX];
  if (file->cycle == -1)
    snprintf(master_file_name, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(master_file_name, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, file->cycle);
  int driver = DB_HDF5;
  DBfile* master = DBCreate(master_file_name, DB_CLOBBER, DB_LOCAL, "Master file", driver);

  // Stick in cycle/time information if needed.
  DBoptlist* optlist = DBMakeOptlist(2);
  if (file->cycle >= 0)
    DBAddOption(optlist, DBOPT_CYCLE, &file->cycle);
  if (file->time != -FLT_MAX)
  {
    double t = (double)file->time;
    DBAddOption(optlist, DBOPT_DTIME, &t);
  }

  int num_files = file->num_files;
  int num_chunks = file->nproc / num_files;

  // Meshes.
  for (int i = 0; i < file->multimeshes->size; ++i)
  {
    multimesh_t* mesh = file->multimeshes->data[i];

    // Mesh.
    char* mesh_names[file->num_files*num_chunks];
    int mesh_types[file->num_files*num_chunks];
    for (int j = 0; j < file->num_files; ++j)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        char mesh_name[FILENAME_MAX];
        mesh_types[num_chunks*j+c] = mesh->type;
        if (file->cycle == -1)
          snprintf(mesh_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, mesh->name);
        else
          snprintf(mesh_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->cycle, c, mesh->name);
        mesh_names[num_chunks*j+c] = string_dup(mesh_name);
      }
    }

    // Write the multimesh.
    int stat = DBPutMultimesh(master, mesh->name, file->num_files*num_chunks, 
                              mesh_names, mesh_types, optlist);
    if (stat == -1)
      polymec_error("Error writing multi-mesh to Silo master file %s.", master_file_name);

    // Clean up.
    for (int j = 0; j < num_files*num_chunks; ++j)
      polymec_free(mesh_names[j]);
  }

  // Variables.
  for (int i = 0; i < file->multivars->size; ++i)
  {
    multivar_t* var = file->multivars->data[i];

    // Fields.
    char* var_names[file->num_files*num_chunks];
    int var_types[num_files*num_chunks];
    for (int j = 0; j < file->num_files; ++j)
    {
      for (int c = 0; c < num_chunks; ++c)
      {
        char var_name[FILENAME_MAX];
        if (file->cycle == -1)
          snprintf(var_name, FILENAME_MAX, "%d/%s.silo:/domain_%d/%s", j, file->prefix, c, var->name);
        else
          snprintf(var_name, FILENAME_MAX, "%d/%s-%d.silo:/domain_%d/%s", j, file->prefix, file->cycle, c, var->name);
        var_names[num_chunks*j+c] = string_dup(var_name);
        var_types[num_chunks*j+c] = var->type;
      }
    }

    // Write the multivariable data.
    DBPutMultivar(master, var->name, num_files*num_chunks, var_names, var_types, optlist);

    // Clean up.
    for (int j = 0; j < num_files*num_chunks; ++j)
      polymec_free(var_names[j]);
  }

  DBFreeOptlist(optlist);
  DBClose(master);
}
#endif

silo_file_t* silo_file_new(MPI_Comm comm,
                           const char* file_prefix,
                           const char* directory,
                           int num_files,
                           int mpi_tag,
                           int cycle,
                           real_t time)
{
  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));

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
  MPI_Comm_size(comm, &file->nproc);
  MPI_Comm_rank(comm, &file->rank);
  if (num_files == -1)
    file->num_files = file->nproc;
  else
    file->num_files = num_files;
  ASSERT(file->num_files <= file->nproc);
  file->mpi_tag = mpi_tag;

  if (file->nproc > 1)
  {
    // We put the entire data set into a directory named after the 
    // prefix, and every process gets its own subdirectory therein.

    // Create the master directory if we need to.
    if ((strcmp(directory, ".") == 0) && (file->nproc > 1))
      polymec_error("silo_file_new: Multi-process filesets cannot be generated in the current working directory.");
    else if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    if (file->rank == 0)
    {
      create_directory(file->directory, S_IRWXU | S_IRWXG);
      MPI_Barrier(comm);
    }
    else
      MPI_Barrier(comm);

    // Initialize poor man's I/O and figure out group ranks.
    file->baton = PMPIO_Init(file->num_files, PMPIO_WRITE, comm, file->mpi_tag, 
        pmpio_create_file, pmpio_open_file, 
        pmpio_close_file, 0);
    file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
    file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

    // Create a subdirectory for each group.
    char group_dir_name[FILENAME_MAX];
    snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
    if (file->rank_in_group == 0)
    {
      create_directory(group_dir_name, S_IRWXU | S_IRWXG);
      MPI_Barrier(comm);
    }
    else
      MPI_Barrier(comm);

    // Determine a file name and directory name.
    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, cycle);
    char silo_dir_name[FILENAME_MAX];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);

    file->multimeshes = ptr_array_new();
    file->multivars = ptr_array_new();
  }
  else
  {
    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

    int driver = DB_HDF5;
    create_directory(file->directory, S_IRWXU | S_IRWXG);
    file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
    DBSetDir(file->dbfile, "/");
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (cycle == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

  int driver = DB_HDF5;
  create_directory(file->directory, S_IRWXU | S_IRWXG);
  file->dbfile = DBCreate(file->filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
  DBSetDir(file->dbfile, "/");
#endif
  file->mode = DB_CLOBBER;
  file->cycle = cycle;
  file->time = time;

  return file;
}

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int num_files,
                            int mpi_tag,
                            int cycle)
{
  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));
  file->mode = DB_READ;
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
  MPI_Comm_size(comm, &file->nproc);
  MPI_Comm_rank(comm, &file->rank);
  if (num_files == -1)
    file->num_files = file->nproc;
  else
    file->num_files = num_files;
  ASSERT(file->num_files <= file->nproc);
  file->mpi_tag = mpi_tag;

  if (file->nproc > 1)
  {
    // We put the entire data set into a directory named after the 
    // prefix, and every process gets its own subdirectory therein.

    // Create the master directory if we need to.
    if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    if (file->rank == 0)
    {
      DIR* master_dir = opendir(file->directory);
      if (master_dir == NULL)
      {
        polymec_error("silo_file_open: Master directory %s does not exist for file prefix %s.",
            file->directory, file->prefix);
      }
      else
        closedir(master_dir);
      MPI_Barrier(comm);
    }
    else
      MPI_Barrier(comm);

    // Initialize poor man's I/O and figure out group ranks.
    file->baton = PMPIO_Init(file->num_files, PMPIO_READ, comm, file->mpi_tag, 
        pmpio_create_file, pmpio_open_file, 
        pmpio_close_file, 0);
    file->group_rank = PMPIO_GroupRank(file->baton, file->rank);
    file->rank_in_group = PMPIO_RankInGroup(file->baton, file->rank);

    // Make sure a subdirectory exists for each group.
    char group_dir_name[FILENAME_MAX];
    snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
    if (file->rank_in_group == 0)
    {
      DIR* group_dir = opendir(group_dir_name);
      if (group_dir == NULL)
      {
        polymec_error("silo_file_open: Group directory %s does not exist for file prefix %s.",
            group_dir_name, file->prefix);
      }
      else
        closedir(group_dir);
      MPI_Barrier(comm);
    }
    else
      MPI_Barrier(comm);

    // Determine a file name and directory name.
    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", group_dir_name, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", group_dir_name, file->prefix, cycle);
    char silo_dir_name[FILENAME_MAX];
    snprintf(silo_dir_name, FILENAME_MAX, "domain_%d", file->rank_in_group);
    file->dbfile = (DBfile*)PMPIO_WaitForBaton(file->baton, file->filename, silo_dir_name);
    file->multimeshes = ptr_array_new();
    file->multivars = ptr_array_new();
  }
  else
  {
    if (strlen(directory) == 0)
      strncpy(file->directory, ".", FILENAME_MAX);
    else
      strncpy(file->directory, directory, FILENAME_MAX);

    if (cycle == -1)
      snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
    else
      snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

    int driver = DB_HDF5;
    file->dbfile = DBOpen(file->filename, driver, file->mode);
    DBSetDir(file->dbfile, "/");
  }
#else
  if (strlen(directory) == 0)
    strncpy(file->directory, ".", FILENAME_MAX);
  else
    strncpy(file->directory, directory, FILENAME_MAX);

  if (cycle == -1)
    snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(file->filename, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, cycle);

  int driver = DB_HDF5;
  file->dbfile = DBOpen(file->filename, driver, file->mode);
  DBSetDir(file->dbfile, "/");
#endif

  // Get cycle/time information.
  if (DBInqVarExists(file->dbfile, "dtime"))
  {
    double dtime;
    DBReadVar(file->dbfile, "dtime", &dtime);
    file->time = (real_t)dtime;
  }
  if (DBInqVarExists(file->dbfile, "cycle"))
    DBReadVar(file->dbfile, "cycle", &file->cycle);

  return file;
}

void silo_file_close(silo_file_t* file)
{
#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
  {
    // Finish working on this process.
    if (file->mode == DB_CLOBBER)
    {
      // Write multi-block objects to the file if needed.
      write_multivars_to_file(file);
    }

    PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
    PMPIO_Finish(file->baton);

    if (file->mode == DB_CLOBBER)
    {
      // Write the uber-master file containing any multiobjects if need be.
      write_master_file(file);
    }

    ptr_array_free(file->multimeshes);
    ptr_array_free(file->multivars);
  }
  else
    DBClose(file->dbfile);
#else
  // Write the file.
  DBClose(file->dbfile);
#endif

  // Clean up.
  polymec_free(file);
}

static void silo_file_write_tags(silo_file_t* file, tagger_t* tagger, const char* tag_list_name)
{
  ASSERT(file->mode == DB_CLOBBER);

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
    DBPutCompoundarray(file->dbfile, tag_list_name, elem_names->data, elem_lengths->data,
                       elem_names->size, tag_data->data, tag_data->size, DB_INT, 0);
  }

  // Clean up.
  int_array_free(elem_lengths);
  string_array_free(elem_names);
  int_array_free(tag_data);
}

static void silo_file_add_multimesh(silo_file_t* file,
                                    const char* mesh_name, 
                                    int silo_mesh_type)
{
  ASSERT(file->mode == DB_CLOBBER);

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
  {
    multimesh_t* mesh = multimesh_new(mesh_name, silo_mesh_type);
    ptr_array_append_with_dtor(file->multimeshes, mesh, DTOR(multimesh_free));
  }
#endif
}

static void silo_file_add_multivar(silo_file_t* file,
                                   const char* mesh_name, 
                                   const char* field_name,
                                   int silo_var_type)
{
  ASSERT((file->mode == DB_CLOBBER) || (file->mode == DB_APPEND));

#if POLYMEC_HAVE_MPI
  if (file->nproc > 1)
  {
    multivar_t* var = multivar_new(mesh_name, field_name, silo_var_type);
    ptr_array_append_with_dtor(file->multivars, var, DTOR(multivar_free));
  }
#endif
}

void silo_file_write_mesh(silo_file_t* file,
                          const char* mesh_name,
                          mesh_t* mesh)
{
  ASSERT(file->mode == DB_CLOBBER);

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
  DBoptlist* optlist = DBMakeOptlist(10);
  char zonelist_name[FILENAME_MAX];
  snprintf(zonelist_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBAddOption(optlist, DBOPT_PHZONELIST, zonelist_name);

  // Write out the 3D polyhedral mesh.
  int num_cells = mesh->num_cells;
  DBPutUcdmesh(file->dbfile, (char*)mesh_name, 3, coordnames, coords,
               num_nodes, num_cells, 0, 0,
               SILO_FLOAT_TYPE, optlist);

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
  int* cell_face_counts = polymec_malloc(sizeof(int) * num_cells);
  memset(cell_face_counts, 0, sizeof(int) * num_cells);
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  DBPutPHZonelist(file->dbfile, zonelist_name, 
                  num_faces, face_node_counts,
                  mesh->face_node_offsets[num_faces], mesh->face_nodes,
                  ext_faces, num_cells, cell_face_counts,
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
  {
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_write_tags(file, mesh->node_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_write_tags(file, mesh->edge_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_face_tags", mesh_name);
    silo_file_write_tags(file, mesh->face_tags, tag_name);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_write_tags(file, mesh->cell_tags, tag_name);
  }

  // Write out the number of mesh cells to a special variable.
  char num_cells_var[FILENAME_MAX];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  int one = 1;
  DBWrite(file->dbfile, num_cells_var, &mesh->num_cells, &one, 1, DB_INT);
  
  // Clean up.
  DBFreeOptlist(optlist);

  // Add a multi-object entry.
  silo_file_add_multimesh(file, mesh_name, DB_UCDMESH);
}

void silo_file_write_scalar_cell_field(silo_file_t* file,
                                       const char* field_name,
                                       const char* mesh_name,
                                       real_t* field_data)
{
  ASSERT(file->mode == DB_CLOBBER);

  // How many cells does our mesh have?
  char num_cells_var[FILENAME_MAX];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  ASSERT(DBInqVarExists(file->dbfile, num_cells_var));
  int num_cells;
  DBReadVar(file->dbfile, num_cells_var, &num_cells);

  // Feed the field data into the file.
  DBPutUcdvar1(file->dbfile, field_name, mesh_name, field_data, num_cells, 0, 0, SILO_FLOAT_TYPE, DB_ZONECENT, NULL);

  // Add a multi-object entry.
  silo_file_add_multivar(file, mesh_name, field_name, DB_UCDVAR);
}

void silo_file_write_cell_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components)
{
  for (int c = 0; c < num_components; ++c)
  {
    silo_file_write_scalar_cell_field(file, field_component_names[c], 
                                      mesh_name, field_data);
  }
}

void silo_file_write_point_mesh(silo_file_t* file,
                                const char* point_mesh_name,
                                point_t* points,
                                int num_points)
{
  ASSERT(file->mode == DB_CLOBBER);

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
  DBPutPointmesh(file->dbfile, (char*)point_mesh_name, 3, coords, num_points, SILO_FLOAT_TYPE, NULL); 
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Write out the number of points to a special variable.
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", point_mesh_name);
  int one = 1;
  DBWrite(file->dbfile, num_points_var, &num_points, &one, 1, DB_INT);
  
  // Add a multi-object entry.
  silo_file_add_multimesh(file, point_mesh_name, DB_POINTMESH);
}

void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* point_mesh_name,
                                        real_t* field_data)
{
  ASSERT(file->mode == DB_CLOBBER);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", point_mesh_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  // Write the point mesh.
  DBPutPointvar1(file->dbfile, field_name, point_mesh_name, field_data, num_points, SILO_FLOAT_TYPE, NULL);

  // Add a multi-object entry.
  silo_file_add_multivar(file, point_mesh_name, field_name, DB_POINTVAR);
}

void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* point_mesh_name,
                                 real_t* field_data,
                                 int num_components)
{
  for (int c = 0; c < num_components; ++c)
  {
    silo_file_write_scalar_point_field(file, field_component_names[c], 
                                       point_mesh_name, field_data);
  }
}

