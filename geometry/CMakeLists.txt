add_library(polymec_geometry voronoi.c planar_voronoi.c 
            plane.c cylinder.c sphere.c intersection.c cubic_lattice.c
            create_cubic_lattice_mesh.c prob_cvt_gen.c 
            interpreter_register_geometry_functions.c)
            #create_cyl_mesh.c create_sphere_mesh.c)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)
