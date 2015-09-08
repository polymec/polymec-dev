#include(add_polymec_test)

#add_mpi_polymec_test(test_partition_point_cloud_with_neighbors test_partition_point_cloud_with_neighbors.c ../../model/tests/create_simple_pairing.c 1 2 3 4)
#add_mpi_polymec_test(test_shepard_shape_function test_shepard_shape_function.c 1 2 3 4)
#add_mpi_polymec_test(test_mls_shape_function test_mls_shape_function.c 1 2 3 4)
#add_polymec_test(test_gmls_functional test_gmls_functional.c poisson_gmls_functional.c make_mlpg_lattice.c)
#add_polymec_test(test_gmls_matrix test_gmls_matrix.c poisson_gmls_functional.c elastic_gmls_functional.c make_mlpg_lattice.c)

