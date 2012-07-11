# Find the PARMETIS includes and libraries
#
# ParMETIS is an MPI-based parallel library that implements a variety of algorithms for 
# partitioning unstructured graphs, meshes, and for computing fill-reducing orderings of 
# sparse matrices. It can be found at:
# 	http://www-users.cs.umn.edu/~karypis/metis/parmetis/index.html
#
# PARMETIS_INCLUDE_DIR - where to find parmetis.h
# PARMETIS_LIBRARIES   - List of fully qualified libraries to link against.
# PARMETIS_FOUND       - Do not attempt to use if "no" or undefined.

# Try to find the ptest executable.
find_program( parmetis_ptest_exe 
    NAMES ptest
    PATH_SUFFIXES bin Bin 
    DOC "Parmetis ptest program." )
mark_as_advanced( parmetis_ptest_exe )

# Intuit directories for Parmetis from ptest.
if (parmetis_ptest_exe)
  string(REPLACE "/bin/ptest" "/include" possible_inc_dir ${parmetis_ptest_exe})
  string(REPLACE "/bin/ptest" "/lib" possible_lib_dir ${parmetis_ptest_exe})
endif ()

FIND_PATH(PARMETIS_INCLUDE_DIR parmetis.h
  /usr/local/include
  /usr/include
  "${possible_inc_dir}"
)

FIND_LIBRARY(PARMETIS_LIBRARY parmetis
  /usr/local/lib
  /usr/lib
  "${possible_lib_dir}"
)

FIND_LIBRARY(METIS_LIBRARY metis
  /usr/local/lib
  /usr/lib
  "${possible_lib_dir}"
)

if (PARMETIS_INCLUDE_DIR)
  if (PARMETIS_LIBRARY)
    SET( PARMETIS_LIBRARIES ${PARMETIS_LIBRARY} ${METIS_LIBRARY})
    SET( PARMETIS_FOUND TRUE)
    message("-- Found Parmetis: ${PARMETIS_LIBRARIES}")
  else()
    SET( PARMETIS_FOUND FALSE)
    message("-- Parmetis not found.")
  endif()
else()
  SET( PARMETIS_FOUND FALSE)
  message("-- Parmetis not found.")
endif()
