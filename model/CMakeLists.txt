add_library(polymec_model model.c boundary_cell_map.c periodic_bc.c stencil.c
            neighbor_pairing.c polynomial_fit.c helmholtz_solver.c conduction_solver.c
            interpreter_register_model_functions.c)

set(POLYMEC_LIBRARIES polymec_model;${POLYMEC_LIBRARIES} PARENT_SCOPE)
set(POLYMEC_LIBRARIES polymec_model;${POLYMEC_LIBRARIES})
add_subdirectory(tests)

# Install library.
install(TARGETS polymec_model DESTINATION lib)

# Install headers, including generated files.
file(GLOB HEADERS *.h)
install(FILES ${HEADERS} DESTINATION include/polymec/model)
