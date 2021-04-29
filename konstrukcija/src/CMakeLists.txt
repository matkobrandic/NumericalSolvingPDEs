add_executable("konstrukcija" yasp_grid.cc)
target_link_dune_default_libraries("konstrukcija")

add_executable("alugrid" alugrid_gmsh.cc)
target_link_dune_default_libraries("konstrukcija")

add_executable("alu_dgf" alugrid_dgf.cc)
target_link_dune_default_libraries("alu_dgf")
