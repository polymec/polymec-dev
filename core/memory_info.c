// Copyright (c) 2012-2015, Jeffrey N. Johnson
// All rights reserved.
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/memory_info.h"

// This stuff is not, ehm, "portable" in any sense of the word. So we do what
// we can.
#ifdef LINUX
void get_memory_info_linux(memory_info_t* info)
{
  info->total_virtual_memory = 0;
  info->virtual_memory_used = 0;
  info->total_physical_memory = 0;
  info->physical_memory_used = 0;
  info->process_virtual_size = 0;
  info->process_resident_size = 0;
  info->process_peak_resident_size = 0;
}
#endif

#ifdef APPLE
#include <mach/mach.h>
#include <mach/vm_statistics.h>
#include <mach/mach_types.h>
#include <mach/mach_init.h>
#include <mach/mach_host.h>
#include <sys/param.h>
#include <sys/mount.h>
#include <sys/sysctl.h>
#include <sys/types.h>

void get_memory_info_apple(memory_info_t* info)
{
  // Available virtual memory.
  struct statfs stats;
  int status = statfs("/", &stats);
  if (status == 0)
  {
    uint64_t total_vm = (uint64_t)stats.f_bsize * stats.f_bfree;
    info->total_virtual_memory = total_vm / 1024;
  }
  else
    info->total_virtual_memory = 0;

  // Used virtual memory.
  struct xsw_usage vmusage;
  info->virtual_memory_used = sizeof(vmusage);
  status = sysctlbyname("vm.swapusage", &vmusage, &info->virtual_memory_used, NULL, 0);
  if (status != 0)
    info->virtual_memory_used = 0;

  // Available physical memory.
  int mib[2];
  int64_t physical_memory;
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  info->total_virtual_memory = sizeof(int64_t);
  sysctl(mib, 2, &physical_memory, &info->total_virtual_memory, NULL, 0);
  info->total_virtual_memory /= 1024;

  // Used physical memory.
  vm_size_t page_size;
  mach_port_t mach_port = mach_host_self();
  vm_statistics64_data_t vm_stats;
  mach_msg_type_number_t count = sizeof(vm_stats) / sizeof(natural_t);
  if ((host_page_size(mach_port, &page_size) == KERN_SUCCESS) &&
      (host_statistics64(mach_port, HOST_VM_INFO,
                         (host_info64_t)&vm_stats, &count) == KERN_SUCCESS))
  {
    long long ram_free = (int64_t)vm_stats.free_count * (int64_t)page_size;
    long long ram_used = ((int64_t)vm_stats.active_count +
                          (int64_t)vm_stats.inactive_count +
                          (int64_t)vm_stats.wire_count) *  (int64_t)page_size;
    info->physical_memory_used = (size_t)(ram_used / 1024);
    info->physical_memory_free = (size_t)(ram_free / 1024);
  }
  else
  {
    info->physical_memory_used = 0;
    info->physical_memory_free = 0;
  }

  // Process memory usage.
  struct task_basic_info t_info;
  mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(),
                TASK_BASIC_INFO, (task_info_t)&t_info, 
                &t_info_count) != KERN_SUCCESS)
  {
    info->process_virtual_size = 0;
    info->process_resident_size = 0;
  }
  else
  {
    info->process_virtual_size = t_info.virtual_size / 1024;
    info->process_resident_size = t_info.resident_size / 1024;
  }

  // Peak resident set size.
  struct rusage rusage;
	getrusage(RUSAGE_SELF, &rusage);
  info->process_peak_resident_size = rusage.ru_maxrss / 1024;

}
#endif

void get_memory_info(memory_info_t* info)
{
#ifdef LINUX
  get_memory_info_linux(info);
#endif
#ifdef APPLE
  get_memory_info_apple(info);
#endif
}

