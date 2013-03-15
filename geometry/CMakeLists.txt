add_library(polymec_geometry plane.c cylinder.c sphere.c intersection.c 
            union.c difference.c scaled.c giftwrap_hull.c polygon2.c polygon.c
            hexahedron.c block.c create_block_mesh.c ball_mapping.c cylinder_mapping.c
            cubic_lattice.c create_cubic_lattice_mesh.c cvt_gen_dist.c
            prob_cvt_gen_dist.c force_cvt_gen_dist.c 
            create_unbounded_voronoi_mesh.c create_bounded_voronoi_mesh.c 
            crop_voronoi_mesh.c voronoi_tessellator.cpp 
            create_manifold.c generate_random_points.c
            delaunay_triangulator.cpp 
            interpreter_register_geometry_functions.c)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS} PARENT_SCOPE)
set(POLYMEC_LIBS polymec_geometry;${POLYMEC_LIBS})

add_subdirectory(tests)
