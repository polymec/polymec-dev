// Copyright (c) 2012-2019, Jeffrey N. Johnson
// All rights reserved.
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "core/memory_info.h"

// This stuff is not, ehm, "portable" in any sense of the word. So we do what
// we can. This stuff is largely taken from the discussion at
// stackoverflow.com/questions/63166/how-to-determine-cpu-and-memory-consumption-from-inside-a-process.
#ifdef LINUX
#include <sys/types.h>
#include <sys/sysinfo.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include "core/logging.h"

static void get_memory_info_linux(memory_info_t* info)
{
  // Available virtual memory.
  struct sysinfo mem_info;
  sysinfo (&mem_info);
  long long vmem = mem_info.totalram;
  vmem += mem_info.totalswap;
  vmem *= mem_info.mem_unit;
  info->total_virtual_memory = vmem / 1024;

  // Used virtual memory.
  long vmem_used = mem_info.totalram - mem_info.freeram;
  vmem_used += mem_info.totalswap - mem_info.freeswap;
  vmem_used *= mem_info.mem_unit;
  info->virtual_memory_used = vmem / 1024;

  // Available physical memory.
  long long total_ram = mem_info.totalram;
  total_ram *= mem_info.mem_unit;
  info->total_physical_memory = total_ram / 1024;

  // Used physical memory.
  long long ram_used = mem_info.totalram - mem_info.freeram;
  ram_used *= mem_info.mem_unit;
  info->physical_memory_used = ram_used / 1024;

  // Free physical memory.
  long long free_ram = mem_info.freeram * mem_info.mem_unit;
  info->physical_memory_free = free_ram / 1024;

  // Process memory usage.
  info->process_virtual_size = 0;
  info->process_resident_size = 0;
  info->process_peak_resident_size = 0;
  FILE* proc = fopen("/proc/self/status", "r");
  static bool first_time = true;
  if (proc != NULL)
  {
    while (!feof(proc))
    {
      char proc_line[80];
      char* proc_ptr = fgets(proc_line, 80, proc);
      if (proc_ptr != NULL)
      {
        // Read off the values of interest (they're already in kB).
        char* p = strstr(proc_line, "VmHWM:");
        if (p != NULL)
        {
          long long unsigned int hwm = 0;
          sscanf(p, "VmHWM:"" %llu", &hwm);
          info->process_peak_resident_size = hwm;
        }
        p = strstr(proc_line, "VmSize:");
        if (p != NULL)
        {
          long long unsigned int size = 0;
          sscanf(p, "VmSize:"" %llu", &size);
          info->process_virtual_size = size;
        }
        p = strstr(proc_line, "VmRSS:");
        if (p != NULL)
        {
          long long unsigned int rss = 0;
          sscanf(p, "VmRSS:"" %llu", &rss);
          info->process_resident_size = rss;
        }
      }
    }
    fclose(proc);
  }
  else if (first_time)
  {
    log_debug("get_memory_info: /proc/self/status is not available.");
    log_debug("get_memory_info: Process memory information will be unreliable.");
  }
  first_time = true;
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

static void get_memory_info_apple(memory_info_t* info)
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
  size_t physical_memory_len = sizeof(int64_t);
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  sysctl(mib, 2, &info->total_physical_memory, &physical_memory_len, NULL, 0);
  info->total_physical_memory /= 1024;

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

