add_executable("integracija" integracija.cc)
target_link_dune_default_libraries("integracija")

add_executable("boundary" bdry_integration.cc)
target_link_dune_default_libraries("boundary")

add_executable("visu" visu1.cc elementdata.hh)
target_link_dune_default_libraries("visu")
