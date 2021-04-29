add_executable("fem" fem.cc)
target_link_dune_default_libraries("fem")

add_custom_command(
        TARGET fem POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/cube.msh
                ${CMAKE_CURRENT_BINARY_DIR}/)
