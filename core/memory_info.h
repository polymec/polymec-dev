// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef POLYMEC_MEMORY_INFO_H
#define POLYMEC_MEMORY_INFO_H

#include <stdlib.h>

/// \addtogroup core core
///@{

/// \struct memory_info
/// This type holds various diagnostics related to memory availability and
/// usage for Linux/UNIX and Mac platforms.
typedef struct memory_info_t
{
  // These fields are all expressed in kB and named so as to be as obvious
  // as possible. :-)
  size_t total_virtual_memory;
  size_t virtual_memory_used;

  size_t total_physical_memory;
  size_t physical_memory_free;
  size_t physical_memory_used;

  size_t process_virtual_size;
  size_t process_resident_size;
  size_t process_peak_resident_size;
} memory_info_t;

/// This function populates the given memory_info struct with data from
/// the system at the time it is called. Fields that can't be filled are
/// set to 0.
/// \relates memory_info
void get_memory_info(memory_info_t* info);

///@}

#endif

