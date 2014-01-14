add_library(polymec_model model.c boundary_cell_map.c 
            periodic_bc.c polynomial_fit.c repartition.c
            interpreter_register_model_functions)

set(POLYMEC_LIBS polymec_model;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_model;${POLYMEC_LIBS})
add_subdirectory(tests)

# Install library.
install(TARGETS polymec_model DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/model)
