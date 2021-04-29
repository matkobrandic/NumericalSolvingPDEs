dune_symlink_to_source_files(FILES "pinch-2D-simplex.msh")
add_executable("tromost" tromost.cc)
target_link_dune_default_libraries("tromost")
