add_library(polymec_geometry voronoi.c planar_voronoi.c 
            plane.c cylinder.c sphere.c intersection.c union.c 
            difference.c giftwrap_hull.c 
            cubic_lattice.c create_cubic_lattice_mesh.c prob_cvt_gen.c 
            interpreter_register_geometry_functions.c)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)
