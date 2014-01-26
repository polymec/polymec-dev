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
#include "model/exchanger.h"

void test_exchanger_new(void** state)
{
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  exchanger_free(exchanger);
}

void test_exchanger_construct(void** state)
{
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  int send_indices[5] = {0, 1, 2, 3, 4};
  exchanger_set_send(exchanger, 1, 5, send_indices, true);
  int receive_indices[5] = {5, 6, 7, 8, 9};
  exchanger_set_receive(exchanger, 1, 5, receive_indices, true);
  exchanger_free(exchanger);
}

void test_exchanger_construct_and_delete(void** state)
{
  exchanger_t* exchanger = exchanger_new(MPI_COMM_WORLD);
  int send_indices[5] = {0, 1, 2, 3, 4};
  exchanger_set_send(exchanger, 1, 5, send_indices, true);
  int receive_indices[5] = {5, 6, 7, 8, 9};
  exchanger_set_receive(exchanger, 1, 5, receive_indices, true);

  exchanger_delete_send(exchanger, 1);
  exchanger_delete_receive(exchanger, 1);
  exchanger_free(exchanger);
}

int main(int argc, char* argv[]) 
{
  polymec_init(argc, argv);
  const UnitTest tests[] = 
  {
    unit_test(test_exchanger_new),
    unit_test(test_exchanger_construct)
  };
  return run_tests(tests);
}
