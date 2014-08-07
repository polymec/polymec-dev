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

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <string.h>
#include "cmockery.h"
#include "core/exchanger.h"

void test_exchanger_new(void** state)
{
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  exchanger_free(exchanger);
}

void test_exchanger_construct(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);
  }
  exchanger_free(exchanger);
}

void test_exchanger_construct_and_delete(void** state)
{
  MPI_Comm comm = MPI_COMM_WORLD;
  int nproc, rank;
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &rank);
  int N = 100*nproc;
  exchanger_t* exchanger = exchanger_new(comm);
  if (nproc > 1)
  {
    int send_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_send(exchanger, (rank+1) % nproc, send_indices, N/nproc, true);
    int receive_indices[N/nproc];
    for (int i = 0; i < N/nproc; ++i)
      send_indices[i] = i;
    exchanger_set_receive(exchanger, (rank+nproc-1) % nproc, receive_indices, N/nproc, true);

    exchanger_delete_send(exchanger, 1);
    exchanger_delete_receive(exchanger, 1);
  }
  exchanger_free(exchanger);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exchanger_new),
    unit_test(test_exchanger_construct),
    unit_test(test_exchanger_construct_and_delete)
  };
  return run_tests(tests);
}
