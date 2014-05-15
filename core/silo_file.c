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

struct silo_file_t 
{
  MPI_Comm comm;
  char prefix[FILENAME_MAX], dir_name[FILENAME_MAX],
       filename[FILENAME_MAX], master_dir_name[FILENAME_MAX],
       group_dir_name[FILENAME_MAX];
  int num_files, mpi_tag, nproc, rank, group_rank, rank_in_group;
  int cycle;
  real_t time;
  PMPIO_baton_t* baton;
  DBfile* dbfile;
};

silo_file_t* silo_file_open(MPI_Comm comm,
                            const char* file_prefix,
                            const char* directory,
                            int num_files,
                            int mpi_tag)
{
  silo_file_t* file = malloc(sizeof(silo_file_t));
  file->comm = comm;
  file->num_files = 1;
  file->mpi_tag = mpi_tag;
  file->nproc = 1;
  file->rank = 0;
  file->group_rank = 0;
  file->rank_in_group = 0;
  file->baton = NULL;
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
#else
  if (strlen(directory) == 0)
    strncpy(file->dir_name, ".", FILENAME_MAX);
  else
    strncpy(file->dir_name, directory, FILENAME_MAX);

  snprintf(file->filename, FILENAME_MAX, "%s/%s.silo", file->dir_name, file->prefix);

  int driver = DB_HDF5;
  file->file = DBCreate(filename, 0, DB_LOCAL, 0, driver);
  DBSetDir(file->file, "/");
#endif

  return file;
}

void silo_file_close(silo_file_t* file)
{
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
}

void silo_file_add_point_mesh(silo_file_t* file,
                              const char* point_mesh_name,
                              point_t* points,
                              int num_points,
                              string_ptr_unordered_map_t* fields)
{
}

