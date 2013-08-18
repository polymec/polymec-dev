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
