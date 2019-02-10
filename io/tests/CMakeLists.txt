include(add_polymec_test)

# These functions create the tests for our io library.
function(add_polymec_io_test exe)
  add_polymec_test_with_libs(${exe} "polymec_io;polymec_model;polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_io_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_io;polymec_model;polymec_solvers;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})

  # Silo is brittle enough with parallel I/O that we don't want more than
  # one of these tests running at a time.
  # NOTE: Here we use some special knowledge from within 
  # NOTE: add_mpi_polymec_test_with_libs.
  foreach (arg ${ARGN})
    if (arg MATCHES ".c")
      list(APPEND sources ${arg})
    else()
      list(APPEND procs ${arg})
    endif()
  endforeach()
  foreach (proc ${procs})
    set_tests_properties(${exe}_${proc}_proc PROPERTIES RUN_SERIAL ON)
  endforeach()
endfunction()

add_mpi_polymec_io_test(test_partition_polymesh_io test_partition_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_partition_point_cloud_io test_partition_point_cloud_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_repartition_polymesh_io test_repartition_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_repartition_point_cloud_io test_repartition_point_cloud_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_create_uniform_polymesh_io test_create_uniform_polymesh_io.c 1 2 4)
add_mpi_polymec_io_test(test_create_rectilinear_polymesh_io test_create_rectilinear_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_crop_polymesh_io test_crop_polymesh_io.c 1 2 4)
add_mpi_polymec_io_test(test_star_stencil_io test_star_stencil_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_silo_file_unimesh_methods test_silo_file_unimesh_methods.c 1 2 4)
add_polymec_io_test(test_planar_polymesh_io test_planar_polymesh_io.c)
add_mpi_polymec_io_test(test_colmesh_io test_colmesh_io.c 1 2 3 4)

