add_executable("test-prog" ug_gmsh.cc)
target_link_dune_default_libraries("test-prog")
dune_symlink_to_source_files(FILES airfoil_exterior-2D-simplex.msh)
