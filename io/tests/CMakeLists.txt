include(add_polymec_test)

# These functions create the tests for our io library.
function(add_polymec_io_test exe)
  add_polymec_test_with_libs(${exe} "polymec_io;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

function(add_mpi_polymec_io_test exe)
  add_mpi_polymec_test_with_libs(${exe} "polymec_io;polymec_geometry;polymec_core;${POLYMEC_BASE_LIBRARIES}" ${ARGN})
endfunction()

add_mpi_polymec_io_test(test_partition_polymesh_io test_partition_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_partition_point_cloud_io test_partition_point_cloud_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_repartition_polymesh_io test_repartition_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_repartition_point_cloud_io test_repartition_point_cloud_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_create_uniform_polymesh_io test_create_uniform_polymesh_io.c 1 2 4)
add_mpi_polymec_io_test(test_create_rectilinear_polymesh_io test_create_rectilinear_polymesh_io.c 1 2 3 4)
add_mpi_polymec_io_test(test_crop_polymesh_io test_crop_polymesh_io.c 1 2 4)


