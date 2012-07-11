if (NOT USE_MPI)
  add_library(mpi_stubs mpi_stubs.c)
  set (BRICK_LIBS ${BRICK_LIBS} mpi_stubs)
endif (NOT USE_MPI)

