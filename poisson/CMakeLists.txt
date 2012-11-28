# Executable.
add_executable(poisson poisson.c laplacian_op.c poisson_model.c)
set_target_properties(poisson PROPERTIES OUTPUT_NAME poisson)
target_link_libraries(poisson ${POLYMEC_LIBS})
#add_subdirectory(tests)
