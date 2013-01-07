# Executable.
add_executable(cnav cnav.c cnav_model.c cnav_bc.c 
               cnav_implicit_model.c cnav_semi_implicit_model.c
               interpreter_register_cnav_functions.c
               register_cnav_benchmarks.c 
               cnav_conduction_solver.c cnav_viscous_solver.c)
set_target_properties(cnav PROPERTIES OUTPUT_NAME cnav)
target_link_libraries(cnav ${POLYMEC_LIBS})
#add_subdirectory(tests)
