// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/polymec.h"
#include "geometry/field_metadata.h"

// Tensor types for field components.
typedef enum
{
  FIELD_UNDEFINED = 0,
  FIELD_VECTOR = 1,
  FIELD_TENSOR2 = 2,
  FIELD_SYMTENSOR2 = 3
} field_comp_type_t;

struct field_metadata_t
{
  // Number of field components.
  int num_comps;
  // Component names.
  char** names; 
  // Component units.
  char** units; 
  // Component conserved quantity flags.
  bool* conserved; 
  // Component extensive quantity flags.
  bool* extensive; 
  // Component types.
  field_comp_type_t* comp_types; 
};

static void field_metadata_free(void* context)
{
  field_metadata_t* md = context;
  for (int c = 0; c < md->num_comps; ++c)
  {
    if (md->names[c] != NULL)
      string_free(md->names[c]);
    if (md->units[c] != NULL)
      string_free(md->units[c]);
  }
  polymec_free(md->names);
  polymec_free(md->units);
  polymec_free(md->conserved);
  polymec_free(md->extensive);
  polymec_free(md->comp_types);
}

field_metadata_t* field_metadata_new(int num_comps)
{
  ASSERT(num_comps > 0);
  field_metadata_t* md = polymec_refcounted_malloc(sizeof(field_metadata_t),
                                                   field_metadata_free);
  md->num_comps = num_comps;
  md->names = polymec_calloc(num_comps, sizeof(char*));
  md->units = polymec_calloc(num_comps, sizeof(char*));
  md->conserved = polymec_calloc(num_comps, sizeof(bool));
  md->extensive = polymec_calloc(num_comps, sizeof(bool));
  md->comp_types = polymec_calloc(num_comps, sizeof(field_comp_type_t));
  return md;
}

field_metadata_t* field_metadata_clone(field_metadata_t* md)
{
  field_metadata_t* clone = field_metadata_new(md->num_comps);
  for (int c = 0; c < md->num_comps; ++c)
  {
    if (md->names[c] != NULL)
      clone->names[c] = string_dup(md->names[c]);
    if (md->units[c] != NULL)
      clone->units[c] = string_dup(md->units[c]);
  }
  memcpy(clone->conserved, md->conserved, sizeof(bool) * md->num_comps);
  memcpy(clone->extensive, md->extensive, sizeof(bool) * md->num_comps);
  memcpy(clone->comp_types, md->comp_types, sizeof(field_comp_type_t) * md->num_comps);
  return clone;
}

int field_metadata_num_components(field_metadata_t* md)
{
  return md->num_comps;
}

const char* field_metadata_name(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  return md->names[component];
}

void field_metadata_set_name(field_metadata_t* md, 
                             int component,
                             const char* name)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  if (md->names[component] != NULL)
    string_free(md->names[component]);
  if (name != NULL)
    md->names[component] = string_dup(name);
  else
    md->names[component] = NULL;
}

const char* field_metadata_units(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  return md->units[component];
}

void field_metadata_set_units(field_metadata_t* md, 
                              int component, 
                              const char* units)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  if (md->units[component] != NULL)
    string_free(md->units[component]);
  if (units != NULL)
    md->units[component] = string_dup(units);
  else
    md->units[component] = NULL;
}

bool field_metadata_conserved(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  return md->conserved[component];
}

void field_metadata_set_conserved(field_metadata_t* md, 
                                  int component,
                                  bool conserved)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  md->conserved[component] = conserved;
}

bool field_metadata_extensive(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  return md->extensive[component];
}

void field_metadata_set_extensive(field_metadata_t* md, 
                                  int component,
                                  bool extensive)
{
  ASSERT(component >= 0);
  ASSERT(component < md->num_comps);
  md->extensive[component] = extensive;
}

void field_metadata_set_vector(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component <= (md->num_comps - 3));
  md->comp_types[component] = FIELD_VECTOR;
}

bool field_metadata_has_vectors(field_metadata_t* md)
{
  for (int c = 0; c < md->num_comps; ++c)
  {
    if (md->comp_types[c] == FIELD_VECTOR)
      return true;
  }
  return false;
}

bool field_metadata_next_vector(field_metadata_t* md,
                                int* pos,
                                int* comp)
{
  if (*pos >= md->num_comps)
    return false;
  while ((*pos < md->num_comps) && (md->comp_types[*pos] != FIELD_VECTOR))
    ++(*pos);
  if (*pos < md->num_comps) 
  {
    ASSERT(md->comp_types[*pos] == FIELD_VECTOR);
    *comp = *pos;
    return true;
  }
  else
    return false;
}

void field_metadata_set_tensor2(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component <= (md->num_comps - 9));
  md->comp_types[component] = FIELD_TENSOR2;
}

bool field_metadata_has_tensor2s(field_metadata_t* md)
{
  for (int c = 0; c < md->num_comps; ++c)
  {
    if (md->comp_types[c] == FIELD_TENSOR2)
      return true;
  }
  return false;
}

bool field_metadata_next_tensor2(field_metadata_t* md,
                                 int* pos,
                                 int* comp)
{
  if (*pos >= md->num_comps)
    return false;
  while ((*pos < md->num_comps) && (md->comp_types[*pos] != FIELD_TENSOR2))
    ++(*pos);
  if (*pos < md->num_comps) 
  {
    ASSERT(md->comp_types[*pos] == FIELD_TENSOR2);
    *comp = *pos;
    return true;
  }
  else
    return false;
}

void field_metadata_set_symtensor2(field_metadata_t* md, int component)
{
  ASSERT(component >= 0);
  ASSERT(component <= (md->num_comps - 6));
  md->comp_types[component] = FIELD_SYMTENSOR2;
}

bool field_metadata_has_symtensor2s(field_metadata_t* md)
{
  for (int c = 0; c < md->num_comps; ++c)
  {
    if (md->comp_types[c] == FIELD_SYMTENSOR2)
      return true;
  }
  return false;
}

bool field_metadata_next_symtensor2(field_metadata_t* md,
                                    int* pos,
                                    int* comp)
{
  if (*pos >= md->num_comps)
    return false;
  while ((*pos < md->num_comps) && (md->comp_types[*pos] != FIELD_SYMTENSOR2))
    ++(*pos);
  if (*pos < md->num_comps) 
  {
    ASSERT(md->comp_types[*pos] == FIELD_TENSOR2);
    *comp = *pos;
    return true;
  }
  else
    return false;
}

