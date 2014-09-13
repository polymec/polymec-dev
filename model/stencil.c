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

#include "model/stencil.h"
#include "model/stencil.h"

stencil_t* stencil_new(const char* name, int num_indices, 
                       int* offsets, int* indices, real_t* weights,
                       exchanger_t* ex)
{
  ASSERT(num_indices > 0);
  stencil_t* s = polymec_malloc(sizeof(stencil_t));
  s->name = string_dup(name);
  s->num_indices = num_indices;
  s->offsets = offsets;
  s->indices = indices;
  s->weights = weights;
  s->ex = ex;
  return s;
}

void stencil_free(stencil_t* stencil)
{
  polymec_free(stencil->name);
  polymec_free(stencil->offsets);
  polymec_free(stencil->indices);
  polymec_free(stencil->weights);
  exchanger_free(stencil->ex);
  polymec_free(stencil);
}

stencil_t* cell_face_stencil_new(mesh_t* mesh)
{
  // Read off the stencil from the mesh.
  int* offsets = polymec_malloc(sizeof(int) * (mesh->num_cells+1));
  memcpy(offsets, mesh->cell_face_offsets, sizeof(int) * (mesh->num_cells+1));
  int* indices = polymec_malloc(sizeof(int) * offsets[mesh->num_cells]);
  real_t* weights = polymec_malloc(sizeof(real_t) * offsets[mesh->num_cells]);
  for (int i = 0; i < offsets[mesh->num_cells]; ++i)
  {
    indices[i] = mesh_face_opp_cell(mesh, mesh->cell_faces[offsets[i]], i);
    weights[i] = 1.0;
  }
  exchanger_t* ex = exchanger_clone(mesh_exchanger(mesh));
  return stencil_new("cell-face stencil", mesh->num_cells, offsets, indices, 
                     weights, ex);
}

stencil_t* cell_edge_stencil_new(mesh_t* mesh)
{
  // FIXME: Todo.
  return NULL;
}

