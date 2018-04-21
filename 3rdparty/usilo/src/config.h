#ifndef CONFIG_H
#define CONFIG_H

// Lean on HDF5's configuration for everything.
#include "H5pubconf.h"

#ifdef H5_HAVE_UNISTD_H
#define HAVE_UNISTD_H H5_HAVE_UNISTD_H
#endif

#ifdef H5_HAVE_STDLIB_H 
#define HAVE_STDLIB_H H5_HAVE_STDLIB_H
#endif

#ifdef H5_HAVE_SYS_TYPES_H 
#define HAVE_SYS_TYPES_H H5_HAVE_SYS_TYPES_H
#endif

#ifdef H5_HAVE_SYS_STAT_H 
#define HAVE_SYS_STAT_H H5_HAVE_SYS_STAT_H
#endif

#define HAVE_SYS_FCNTL_H 1

#ifdef H5_HAVE_FCNTL_H 
#define HAVE_FCNTL_H H5_HAVE_FCNTL_H
#endif

#ifdef H5_HAVE_STRING_H 
#define HAVE_STRING_H H5_HAVE_STRING_H
#endif

#ifdef H5_HAVE_STRINGS_H 
#define HAVE_STRINGS_H H5_HAVE_STRINGS_H
#endif

#define SIZEOF_OFF64_T 4 // FIXME?

#define HAVE_STRERROR 1

#define HAVE_HDF5_H 1
#define HAVE_LIBHDF5 1
#define HAVE_HDF5_DRIVER 1
#define HAVE_LIBZ 1

#define PACKAGE "usilo"
#define PACKAGE_BUGREPORT "--"
#define PACKAGE_NAME "usilo"
#define PACKAGE_STRING "usilo @USILO_VERSION@"
#define PACKAGE_TARNAME "usilo"
#define PACKAGE_VERSION "@USILO_VERSION@"

#endif
