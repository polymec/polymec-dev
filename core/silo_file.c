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
#include "core/array_utils.h"

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
  DBfile* file = DBCreate(filename, DB_CLOBBER, DB_LOCAL, NULL, driver);
  if (strcmp(dir_name, "/") != 0)
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
    if (strcmp(dir_name, "/") != 0)
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

bool silo_file_query(const char* file_prefix,
                     const char* directory,
                     int* num_files,
                     int* num_mpi_processes,
                     int_slist_t* cycles)
{
  // No blank strings allowed for queries.
  ASSERT(strlen(file_prefix) > 0);
  ASSERT(strlen(directory) > 0);

  // Rank 0 does all the dirty work.
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0)
  {
    // Inspect the given directory's contents.
    if (!directory_exists(directory))
    {
      log_urgent("Can't query non-existent directory: %s", directory);
      return false;
    }
    string_slist_t* files_in_dir = files_within_directory(directory);

    // Try to find a master file or single data file.
    char data_file[FILENAME_MAX];
    data_file[0] = '\0';
    bool found_cycles = false;
    {
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX];
        snprintf(path, FILENAME_MAX, "%s", file_prefix); 
        if (strstr(node->value, path) && strstr(node->value, "silo"))
        {
          strncpy(data_file, node->value, FILENAME_MAX);
          snprintf(data_file, FILENAME_MAX, "%s/%s", directory, node->value);
          snprintf(path, FILENAME_MAX, "%s-", file_prefix); 
          if (strstr(node->value, path))
            found_cycles = true;
          break;
        }
        node = node->next;
      }
    }
    if (strlen(data_file) == 0)
    {
      log_urgent("silo_file_query: Could not find %s/%s-*.silo.", directory, file_prefix);
      string_slist_free(files_in_dir);
      return false;
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
      if (!DBInqVarExists(file, "POLYMEC_SILO_MASTER_FILE"))
      {
        log_urgent("silo_query_file: invalid Silo master file.");
        DBClose(file);
        return false;
      }

      // How many MPI processes were used to construct the data set?
      int my_num_mpi_procs = -1;
      for (int f = 0; f < toc->nmultimesh; ++f)
      {
        DBmultimesh* multimesh = DBGetMultimesh(file, toc->multimesh_names[f]);
        if (my_num_mpi_procs == -1)
          my_num_mpi_procs = multimesh->nblocks;
        ASSERT(my_num_mpi_procs == multimesh->nblocks);
      }
      *num_mpi_processes = my_num_mpi_procs;

      // How many files are in the data set?
      if (DBInqVarExists(file, "num_files"))
        DBReadVar(file, "num_files", num_files);
      else
      {
        log_urgent("silo_file_query: Could not read number of files in set.");
        DBClose(file);
        return false;
      }
    }
    else
    {
      if (!DBInqVarExists(file, "POLYMEC_SILO_FILE"))
      {
        log_urgent("silo_query_file: Invalid Silo file.");
        DBClose(file);
        return false;
      }

      // A single data file can only be written for a serial run.
      *num_files = 1;
      *num_mpi_processes = 1;
    }

    DBClose(file);

    // Search for available cycles.
    if ((cycles != NULL) && found_cycles)
    {
      int_slist_clear(cycles);
      string_slist_node_t* node = files_in_dir->front;
      while (node != NULL)
      {
        char path[FILENAME_MAX];
        snprintf(path, FILENAME_MAX, "%s-", file_prefix); 
        char* p1 = strstr(node->value, path);
        char* p2 = strstr(node->value, ".silo");
        if ((p1 != NULL) && (p2 != NULL))
        {
          char* c = p1 + strlen(path);
          char num[p2-c+1];
          strncpy(num, c, p2-c);
          num[p2-c] = '\0';
          if (string_is_number(num))
            int_slist_append(cycles, atoi(num));
        }
        node = node->next;
      }
    }

    // Clean up.
    string_slist_free(files_in_dir);
  }

  // Now spread the word to other processes.
  int num_cycles = (cycles != NULL) ? cycles->size : 0;
  int data[3] = {*num_files, *num_mpi_processes, num_cycles};
  MPI_Bcast(data, 3, MPI_INT, 0, MPI_COMM_WORLD);
  if (rank != 0)
  {
    *num_files = data[0];
    *num_mpi_processes = data[1];
    num_cycles = data[2];
  }

  if (cycles != NULL)
  {
    int cycles_buffer[num_cycles];
    if (rank == 0)
    {
      // Spit the cycles into the array and sort them.
      int_slist_node_t* node = cycles->front;
      int i = 0;
      while (node != NULL)
      {
        cycles_buffer[i++] = node->value;
        node = node->next;
      }
      ASSERT(i == num_cycles);
      int_qsort(cycles_buffer, num_cycles);
    }
    MPI_Bcast(cycles_buffer, num_cycles, MPI_INT, 0, MPI_COMM_WORLD);

    // Shuffle them back into our linked list.
    int_slist_clear(cycles);
    for (int i = 0; i < num_cycles; ++i)
      int_slist_append(cycles, cycles_buffer[i]);
  }
  return true;
}

struct silo_file_t 
{
  // File data.
  DBfile* dbfile;

  // Metadata.
  char prefix[FILENAME_MAX], directory[FILENAME_MAX], filename[FILENAME_MAX];
  int cycle;
  real_t time;
  int mode; // Open for reading (DB_READ) or writing (DB_CLOBBER)? 
  string_ptr_unordered_map_t* expressions;

#if POLYMEC_HAVE_MPI
  // Stuff for poor man's parallel I/O.
  PMPIO_baton_t* baton;
  MPI_Comm comm;
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;
  ptr_array_t* multimeshes;
  ptr_array_t* multivars;
#endif
};

// Expression struct.
typedef struct
{
  int type;
  char* definition;
} silo_expression_t;

static void write_expressions_to_file(silo_file_t* file, DBfile* dbfile)
{
  ASSERT(file->mode == DB_CLOBBER);
  ASSERT(file->expressions != NULL);

  // Write out the expressions.
  int pos = 0, k = 0;
  char* name;
  void* val;
  char* names[file->expressions->size];
  int types[file->expressions->size];
  char* defs[file->expressions->size];
  while (string_ptr_unordered_map_next(file->expressions, &pos, &name, &val))
  {
    names[k] = name;
    silo_expression_t* exp = val;
    types[k] = exp->type;
    defs[k] = exp->definition;
    ++k;
  }
  DBPutDefvars(dbfile, "expressions", file->expressions->size,
               (const char* const*)names, (const int*)types, 
               (const char* const*)defs, NULL);
}

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
    DBPutMultimesh(file->dbfile, mesh->name, num_chunks, 
                   (char const* const*)mesh_names, mesh_types, optlist);

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
    DBPutMultivar(file->dbfile, var->name, num_chunks, 
                  (char const* const*)var_names, var_types, optlist); 

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

  char master_file_name[FILENAME_MAX];
  if (file->cycle == -1)
    snprintf(master_file_name, FILENAME_MAX, "%s/%s.silo", file->directory, file->prefix);
  else
    snprintf(master_file_name, FILENAME_MAX, "%s/%s-%d.silo", file->directory, file->prefix, file->cycle);
  PMPIO_baton_t* baton = PMPIO_Init(file->num_files, PMPIO_WRITE, file->comm, file->mpi_tag+1, 
                                    pmpio_create_file, pmpio_open_file, 
                                    pmpio_close_file, 0);
  DBfile* master = (DBfile*)PMPIO_WaitForBaton(baton, master_file_name, "/");

  // Write our stamp of approval.
  int one = 1;
  DBWrite(master, "POLYMEC_SILO_MASTER_FILE", &one, &one, 1, DB_INT);

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
                              (char const* const*)mesh_names, mesh_types, 
                              optlist);
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
    DBPutMultivar(master, var->name, num_files*num_chunks, 
                  (char const* const *)var_names, var_types, optlist);

    // Finally, write the number of files to the master file.
    int one = 1;
    DBWrite(master, "num_files", &num_files, &one, 1, DB_INT);

    // Clean up.
    for (int j = 0; j < num_files*num_chunks; ++j)
      polymec_free(var_names[j]);
  }

  // Don't forget to write expressions to the master file.
  write_expressions_to_file(file, master);

  DBFreeOptlist(optlist);

  PMPIO_HandOffBaton(baton, (void*)master);
  PMPIO_Finish(baton);
//  DBClose(master);
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
  file->expressions = string_ptr_unordered_map_new();

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
  file->comm = comm;
  MPI_Comm_size(file->comm, &file->nproc);
  MPI_Comm_rank(file->comm, &file->rank);
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
      if (strcmp(file->directory, ".") != 0)
        create_directory(file->directory, S_IRWXU | S_IRWXG);
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
    char group_dir_name[FILENAME_MAX];
    snprintf(group_dir_name, FILENAME_MAX, "%s/%d", file->directory, file->group_rank);
    if (file->rank_in_group == 0)
    {
      create_directory(group_dir_name, S_IRWXU | S_IRWXG);
      MPI_Barrier(file->comm);
    }
    else
      MPI_Barrier(file->comm);

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
    if (strcmp(file->directory, ".") != 0)
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

  // Write our stamp of approval.
  int one = 1;
  DBWrite(file->dbfile, "POLYMEC_SILO_FILE", &one, &one, 1, DB_INT);

  return file;
}

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int mpi_tag,
                            int cycle, 
                            real_t* time)
{
  silo_file_t* file = polymec_malloc(sizeof(silo_file_t));
  file->mode = DB_READ;
  file->cycle = -1;
  file->time = -FLT_MAX;
  file->expressions = NULL;

  // Strip .silo off of the prefix if it's there.
  {
    char prefix[FILENAME_MAX];
    strncpy(prefix, file_prefix, FILENAME_MAX);
    char* suffix = strstr(prefix, ".silo");
    if (suffix != NULL)
      suffix[0] = '\0';
    strcpy(file->prefix, prefix);
  }

  // Query the dataset for the number of files and MPI processes and cycles.
  int num_files, num_mpi_procs;
  int_slist_t* cycles = int_slist_new();
  if (!silo_file_query(file_prefix, directory, &num_files, &num_mpi_procs, cycles))
    polymec_error("silo_file_open: Invalid file.");

  // For now, we only support reading files that were written with the same 
  // number of processes.
  int nproc = 1;
#if POLYMEC_HAVE_MPI
  MPI_Comm_rank(file->comm, &nproc); 
#endif
  if (nproc != num_mpi_procs)
    polymec_not_implemented("silo_file_open: reading files written with different number of MPI processes is not yet supported.");

  // Check to see whether the requested cycle is available, or whether the 
  // latest one is requested (with -1).
  if (cycle >= 0)
  {
    bool cycle_found = false;
    int_slist_node_t* node = cycles->front;
    while (node != NULL)
    {
      if (node->value == cycle)
      {
        cycle_found = true;
        break;
      }
      else if (node->value > cycle) // cycles are sorted
        break;
      node = node->next;
    }
    if (!cycle_found)
      polymec_error("silo_file_open: Cycle %d was not found for prefix '%s' in directory %s.", cycle, file->prefix, directory);
  }

#if POLYMEC_HAVE_MPI

  // The way these things are defined for a file has to do with how the 
  // file was generated, not how we are currently running.
  file->comm = comm; // ...for lack of a better value. Plus, might be useful.
  MPI_Comm_rank(file->comm, &file->rank); // ...also might be useful.
  file->num_files = num_files; // number of files in the data set.
  file->nproc = num_mpi_procs; // number of MPI procs used to write the thing.
  file->mpi_tag = mpi_tag; // this is fine.

  if (file->nproc > 1)
  {
    // Look for the master directory.
    if (strlen(directory) == 0)
      snprintf(file->directory, FILENAME_MAX, "%s_%dprocs", file->prefix, file->nproc);
    else
      strncpy(file->directory, directory, FILENAME_MAX);
    if (file->rank == 0)
    {
      if (!directory_exists(file->directory))
      {
        polymec_error("silo_file_open: Master directory %s does not exist for file prefix %s.",
            file->directory, file->prefix);
      }
    }

    // Initialize poor man's I/O and figure out group ranks.
    MPI_Barrier(file->comm);
    file->baton = PMPIO_Init(file->num_files, PMPIO_READ, file->comm, file->mpi_tag, 
                             pmpio_create_file, pmpio_open_file, pmpio_close_file, 0);
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
      MPI_Barrier(file->comm);
    }
    else
      MPI_Barrier(file->comm);

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
  else
    file->time = 0.0;
  if (DBInqVarExists(file->dbfile, "cycle"))
    DBReadVar(file->dbfile, "cycle", &file->cycle);
  else
    file->cycle = -1;

  if (time != NULL)
    *time = file->time;

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
      write_expressions_to_file(file, file->dbfile);
    }

    PMPIO_HandOffBaton(file->baton, (void*)file->dbfile);
    PMPIO_Finish(file->baton);

    if (file->mode == DB_CLOBBER)
    {
      // Write the uber-master file containing any multiobjects if need be.
      write_master_file(file);
    }
    MPI_Barrier(file->comm);

    ptr_array_free(file->multimeshes);
    ptr_array_free(file->multivars);
  }
  else
  {
    if (file->mode == DB_CLOBBER)
      write_expressions_to_file(file, file->dbfile);
    DBClose(file->dbfile);
  }
#else
  // Write the file.
  if (file->mode == DB_CLOBBER)
    write_expressions_to_file(file, file->dbfile);
  DBClose(file->dbfile);
#endif

  // Clean up.
  if (file->expressions != NULL)
    string_ptr_unordered_map_free(file->expressions);
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
    DBPutCompoundarray(file->dbfile, tag_list_name, 
                       (char const* const*)elem_names->data, elem_lengths->data,
                       elem_names->size, tag_data->data, tag_data->size, DB_INT, 0);
  }

  // Clean up.
  int_array_free(elem_lengths);
  string_array_free(elem_names);
  int_array_free(tag_data);
}

static void silo_file_read_tags(silo_file_t* file, const char* tag_list_name, tagger_t* tagger)
{
  ASSERT(DBInqVarExists(file->dbfile, tag_list_name));
  DBcompoundarray* var = DBGetCompoundarray(file->dbfile, (char*)tag_list_name);
  int num_tags = var->nelems;
  char** tag_names = var->elemnames;
  int* tag_sizes = var->elemlengths;
  int* array = var->values;
  int size = var->nvalues;
  int j = 0;
  for (int i = 0; i < num_tags; ++i)
  {
    int* tag = tagger_create_tag(tagger, (const char*)tag_names[i], tag_sizes[i]);
    memcpy(tag, &array[j], sizeof(int) * tag_sizes[i]);
    j += tag_sizes[i];
  }
  ASSERT(j == size);
}

static void silo_file_write_exchanger(silo_file_t* file, exchanger_t* ex, const char* exchanger_name)
{
  // Collapse the exchanger into a set of integers.
  // Format is: [nprocs rank send_map receive_map]
  // where send_map and receive_map are sets of integers encoding mappings
  // of processes to sets of indices. Such a map has the format
  // [num_procs [process num_indices i0 i1 ... iN] ... ]
  int_array_t* array = int_array_new();
#if POLYMEC_HAVE_MPI
  int_array_append(array, file->nproc);
  int_array_append(array, file->rank);
#endif

  int pos = 0, proc, *indices, num_indices;
  int_array_append(array, exchanger_num_sends(ex));
  while (exchanger_next_send(ex, &pos, &proc, &indices, &num_indices))
  {
    int_array_append(array, proc);
    int_array_append(array, num_indices);
    for (int i = 0; i < num_indices; ++i)
    int_array_append(array, indices[i]);
  }
  pos = 0;
  int_array_append(array, exchanger_num_receives(ex));
  while (exchanger_next_receive(ex, &pos, &proc, &indices, &num_indices))
  {
    int_array_append(array, proc);
    int_array_append(array, num_indices);
    for (int i = 0; i < num_indices; ++i)
    int_array_append(array, indices[i]);
  }

  // Write the exchanger array to the file.
  int size = array->size;
  DBWrite(file->dbfile, exchanger_name, array->data, &size, 1, DB_INT);
  char exchanger_size_name[FILENAME_MAX];
  snprintf(exchanger_size_name, FILENAME_MAX, "%s_size", exchanger_name);
  int one = 1;
  DBWrite(file->dbfile, exchanger_size_name, &array->size, &one, 1, DB_INT);

  // Clean up.
  int_array_free(array);
}

static void silo_file_read_exchanger(silo_file_t* file, const char* exchanger_name, exchanger_t* ex)
{
  ASSERT(DBInqVarExists(file->dbfile, exchanger_name));
  char exchanger_size_name[FILENAME_MAX];
  snprintf(exchanger_size_name, FILENAME_MAX, "%s_size", exchanger_name);
  ASSERT(DBInqVarExists(file->dbfile, exchanger_name));
  int size;
  DBReadVar(file->dbfile, (char*)exchanger_size_name, &size);
  int* array = polymec_malloc(sizeof(int) * size);
  DBReadVar(file->dbfile, (char*)exchanger_name, array);

  // Expand the exchanger.
  int i = 0;
#if POLYMEC_HAVE_MPI
  ASSERT(array[i++] == file->nproc);
  ASSERT(array[i++] == file->nproc);
#endif
  int num_sends = array[i++];
  for (int j = 0; j < num_sends; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_send(ex, proc, &array[i], num_indices, true);
    array += num_indices;
  }
  int num_receives = array[i++];
  for (int j = 0; j < num_receives; ++j)
  {
    int proc = array[i++];
    int num_indices = array[i++];
    exchanger_set_receive(ex, proc, &array[i], num_indices, true);
    array += num_indices;
  }
  ASSERT(i == size);

  // Clean up.
  polymec_free(array);
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
  DBPutUcdmesh(file->dbfile, (char*)mesh_name, 3, (char const* const*)coordnames, coords,
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
  int* cell_face_counts = polymec_malloc(sizeof(int) * (num_cells + mesh->num_ghost_cells));
  memset(cell_face_counts, 0, sizeof(int) * (num_cells + mesh->num_ghost_cells));
  for (int i = 0; i < num_cells; ++i)
    cell_face_counts[i] = mesh->cell_face_offsets[i+1] - mesh->cell_face_offsets[i];

  // Write the connectivity information.
  DBPutPHZonelist(file->dbfile, zonelist_name, num_faces, face_node_counts,
                  mesh->face_node_offsets[num_faces], mesh->face_nodes,
                  ext_faces, num_cells + mesh->num_ghost_cells, cell_face_counts,
                  mesh->cell_face_offsets[num_cells], mesh->cell_faces,
                  0, 0, num_cells-1, optlist);

  // Partial cleanup.
  polymec_free(face_node_counts);
  polymec_free(ext_faces);
  polymec_free(cell_face_counts);

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

  // Write out exchanger information.
  {
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_exchanger", mesh_name);
    silo_file_write_exchanger(file, mesh_exchanger(mesh), ex_name);
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

mesh_t* silo_file_read_mesh(silo_file_t* file,
                            const char* mesh_name)
{
  ASSERT(file->mode == DB_READ);

  DBucdmesh* ucd_mesh = DBGetUcdmesh(file->dbfile, mesh_name);
  if (ucd_mesh == NULL)
    polymec_error("No mesh named '%s' was found within the Silo file.", mesh_name);
  ASSERT(ucd_mesh->ndims == 3);

  // Also get the polyhedral zone list.
  char phzl_name[FILENAME_MAX];
  snprintf(phzl_name, FILENAME_MAX, "%s_zonelist", mesh_name);
  DBphzonelist* ph_zonelist = DBGetPHZonelist(file->dbfile, phzl_name);
  if (ph_zonelist == NULL)
    polymec_error("Mesh '%s' is not a polymec polyhedral mesh.", mesh_name);

  // Decipher the mesh object.
  int num_cells = ph_zonelist->hi_offset;
  int num_ghost_cells = ph_zonelist->nzones - num_cells;
  int num_faces = ph_zonelist->nfaces;
  int num_nodes = ucd_mesh->nnodes;
#if POLYMEC_HAVE_MPI
  MPI_Comm comm = file->comm;
#else
  MPI_Comm comm = MPI_COMM_WORLD;
#endif
  mesh_t* mesh = mesh_new(comm, num_cells, num_ghost_cells, 
                          num_faces, num_nodes);

  // Set node positions.
  double* x = ucd_mesh->coords[0];
  double* y = ucd_mesh->coords[1];
  double* z = ucd_mesh->coords[2];
  for (int n = 0; n < num_nodes; ++n)
  {
    mesh->nodes[n].x = x[n];
    mesh->nodes[n].y = y[n];
    mesh->nodes[n].z = z[n];
  }

  // Set up cell face counts and face node counts.
  mesh->cell_face_offsets[0] = 0;
  for (int c = 0; c < num_cells; ++c)
    mesh->cell_face_offsets[c+1] = ph_zonelist->facecnt[c];
  mesh->face_node_offsets[0] = 0;
  for (int f = 0; f < num_faces; ++f)
    mesh->face_node_offsets[f+1] = ph_zonelist->nodecnt[f];
  mesh_reserve_connectivity_storage(mesh);

  // Fill in cell faces and face nodes.
  memcpy(mesh->cell_faces, ph_zonelist->facelist, sizeof(int) * mesh->cell_face_offsets[mesh->num_cells]);
  memcpy(mesh->face_nodes, ph_zonelist->nodelist, sizeof(int) * mesh->face_node_offsets[mesh->num_faces]);

  // Finish constructing the mesh.
  mesh_construct_edges(mesh);
  mesh_compute_geometry(mesh);

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->node_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_edge_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->edge_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_face_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->face_tags);
    snprintf(tag_name, FILENAME_MAX, "%s_cell_tags", mesh_name);
    silo_file_read_tags(file, tag_name, mesh->cell_tags);
  }

  // Read in exchanger information.
  {
    char ex_name[FILENAME_MAX];
    snprintf(ex_name, FILENAME_MAX, "%s_exchanger", ex_name);
    silo_file_read_exchanger(file, ex_name, mesh_exchanger(mesh));
  }

  return mesh;
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

real_t* silo_file_read_scalar_cell_field(silo_file_t* file,
                                         const char* field_name,
                                         const char* mesh_name)
{
  ASSERT(file->mode == DB_READ);

  DBucdvar* var = DBGetUcdvar(file->dbfile, (char*)field_name);
  if (var == NULL)
    polymec_error("Field '%s' was not found in the Silo file.", field_name);
  if (var->centering != DB_ZONECENT)
    polymec_error("Field '%s' is not a polymec cell-centered field.", field_name);
  real_t* field = polymec_malloc(sizeof(real_t) * var->nels);
  memcpy(field, var->vals[0], sizeof(real_t) * var->nels);
  return field;
}

void silo_file_write_cell_field(silo_file_t* file,
                                const char** field_component_names,
                                const char* mesh_name,
                                real_t* field_data,
                                int num_components)
{
  ASSERT(file->mode == DB_CLOBBER);

  // How many cells does our mesh have?
  char num_cells_var[FILENAME_MAX];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  ASSERT(DBInqVarExists(file->dbfile, num_cells_var));
  int num_cells;
  DBReadVar(file->dbfile, num_cells_var, &num_cells);
  real_t* comp_data = polymec_malloc(sizeof(real_t) * num_cells); 
  for (int c = 0; c < num_components; ++c)
  {
    for (int i = 0; i < num_cells; ++i)
      comp_data[i] = field_data[num_components*i+c];
    silo_file_write_scalar_cell_field(file, field_component_names[c], 
                                      mesh_name, comp_data);
  }
  polymec_free(comp_data);
}

real_t* silo_file_read_cell_field(silo_file_t* file,
                                  const char** field_component_names,
                                  const char* mesh_name,
                                  int num_components)
{
  ASSERT(file->mode == DB_READ);

  // How many cells does our mesh have?
  char num_cells_var[FILENAME_MAX];
  snprintf(num_cells_var, FILENAME_MAX, "%s_mesh_num_cells", mesh_name);
  ASSERT(DBInqVarExists(file->dbfile, num_cells_var));
  int num_cells;
  DBReadVar(file->dbfile, num_cells_var, &num_cells);
  real_t* field = polymec_malloc(sizeof(real_t) * num_components * num_cells); 
  for (int c = 0; c < num_components; ++c)
  {
    real_t* comp_data = silo_file_read_scalar_cell_field(file, field_component_names[c], mesh_name);
    for (int i = 0; i < num_cells; ++i)
      field[num_components*i+c] = comp_data[i];
    polymec_free(comp_data);
  }
  return field;
}

void silo_file_write_point_cloud(silo_file_t* file,
                                 const char* cloud_name,
                                 point_cloud_t* cloud)
{
  ASSERT(file->mode == DB_CLOBBER);

  // Point coordinates.
  int num_points = cloud->num_points;
  point_t* points = cloud->points;
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
  DBPutPointmesh(file->dbfile, (char*)cloud_name, 3, coords, num_points, SILO_FLOAT_TYPE, NULL); 
  polymec_free(x);
  polymec_free(y);
  polymec_free(z);

  // Write out the number of points to a special variable.
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  int one = 1;
  DBWrite(file->dbfile, num_points_var, &num_points, &one, 1, DB_INT);
  
  // Write out tag information.
  {
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_write_tags(file, cloud->tags, tag_name);
  }

  // Add a multi-object entry.
  silo_file_add_multimesh(file, cloud_name, DB_POINTMESH);
}

point_cloud_t* silo_file_read_point_cloud(silo_file_t* file,
                                          const char* cloud_name)
{
  ASSERT(file->mode == DB_READ);

  // How many points does our cloud have?
  int num_points;
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  DBReadVar(file->dbfile, num_points_var, &num_points);

  DBpointmesh* pm = DBGetPointmesh(file->dbfile, (char*)cloud_name);
  if (pm == NULL)
    polymec_error("Point mesh '%s' was not found in the Silo file.", cloud_name);
  point_t* points = polymec_malloc(sizeof(point_t) * num_points);
  double* x = pm->coords[0];
  double* y = pm->coords[1];
  double* z = pm->coords[2];
  for (int p = 0; p < num_points; ++p)
  {
    points[p].x = x[p];
    points[p].y = y[p];
    points[p].z = z[p];
  }
#if POLYMEC_HAVE_MPI
  point_cloud_t* cloud = point_cloud_from_points(file->comm, points, num_points);
#else
  point_cloud_t* cloud = point_cloud_from_points(MPI_COMM_WORLD, points, num_points);
#endif
  polymec_free(points);

  // Read in tag information.
  {
    char tag_name[FILENAME_MAX];
    snprintf(tag_name, FILENAME_MAX, "%s_node_tags", cloud_name);
    silo_file_read_tags(file, tag_name, cloud->tags);
  }

  return cloud;
}

void silo_file_write_scalar_point_field(silo_file_t* file,
                                        const char* field_name,
                                        const char* cloud_name,
                                        real_t* field_data)
{
  ASSERT(file->mode == DB_CLOBBER);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);

  // Write the field.
  DBPutPointvar1(file->dbfile, field_name, cloud_name, field_data, num_points, SILO_FLOAT_TYPE, NULL);

  // Add a multi-object entry.
  silo_file_add_multivar(file, cloud_name, field_name, DB_POINTVAR);
}

real_t* silo_file_read_scalar_point_field(silo_file_t* file,
                                          const char* field_name,
                                          const char* cloud_name)
{
  ASSERT(file->mode == DB_READ);

  DBmeshvar* var = DBGetPointvar(file->dbfile, (char*)field_name);
  if (var == NULL)
    polymec_error("Field '%s' was not found in the Silo file.", field_name);
  real_t* field = polymec_malloc(sizeof(real_t) * var->nels);
  memcpy(field, var->vals[0], sizeof(real_t) * var->nels);
  return field;
}

void silo_file_write_point_field(silo_file_t* file,
                                 const char** field_component_names,
                                 const char* cloud_name,
                                 real_t* field_data,
                                 int num_components)
{
  ASSERT(file->mode == DB_CLOBBER);

  for (int c = 0; c < num_components; ++c)
  {
    silo_file_write_scalar_point_field(file, field_component_names[c], 
                                       cloud_name, field_data);
  }
}

real_t* silo_file_read_point_field(silo_file_t* file,
                                   const char** field_component_names,
                                   const char* cloud_name,
                                   int num_components)
{
  ASSERT(file->mode == DB_READ);

  // How many points does our mesh have?
  char num_points_var[FILENAME_MAX];
  snprintf(num_points_var, FILENAME_MAX, "%s_num_points", cloud_name);
  ASSERT(DBInqVarExists(file->dbfile, num_points_var));
  int num_points;
  DBReadVar(file->dbfile, num_points_var, &num_points);
  real_t* field = polymec_malloc(sizeof(real_t) * num_components * num_points); 
  for (int c = 0; c < num_components; ++c)
  {
    real_t* comp_data = silo_file_read_scalar_cell_field(file, field_component_names[c], cloud_name);
    for (int i = 0; i < num_points; ++i)
      field[num_components*i+c] = comp_data[i];
    polymec_free(comp_data);
  }
  return field;
}

static void expression_dtor(char* key, void* val)
{
  string_free(key);
  silo_expression_t* exp = val;
  string_free(exp->definition);
  polymec_free(exp);
}

void silo_file_write_scalar_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)

{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_SCALAR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}

void silo_file_write_vector_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)
{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_VECTOR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}

void silo_file_write_tensor_expression(silo_file_t* file,
                                       const char* expression_name,
                                       const char* definition)
{
  ASSERT(file->mode == DB_CLOBBER);
  silo_expression_t* exp = polymec_malloc(sizeof(silo_expression_t));
  exp->type = DB_VARTYPE_TENSOR;
  exp->definition = string_dup(definition);
  string_ptr_unordered_map_insert_with_kv_dtor(file->expressions,
                                               string_dup(expression_name), exp,
                                               expression_dtor);
}
