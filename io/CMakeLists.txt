set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -DPOLYMEC_C_COMPILER=${CMAKE_C_COMPILER} -DPOLYMEC_INCLUDE_DIR=${PROJECT_SOURCE_DIR}")
add_library(polymec_io silo_io.c gnuplot_io.c
            generate_face_node_conn.c generate_cell_node_conn.c)
set(POLYMEC_LIBS polymec_io;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_io;${POLYMEC_LIBS})

add_subdirectory(tests)

