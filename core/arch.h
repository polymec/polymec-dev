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

#ifndef POLYMEC_ARCH_H
#define POLYMEC_ARCH_H

#include <stdio.h>

// This file contains functions that are unavailable on certain target 
// architectures.

#ifdef APPLE

// This is a port of the fmemopen function (available on Linux).
FILE* fmemopen(void *buf, size_t size, const char *mode);

// This is a port of open_memstream (available on Linux).
FILE* open_memstream(char **buf, size_t *len);

#endif 
#endif
