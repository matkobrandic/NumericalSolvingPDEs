add_executable("min-surface" min-surface.cc driver.hh  driver_simple_iter.hh
               l2norm.hh operator.hh  operator_simple_iter.hh bctype.hh
               nlin.ini)
target_link_dune_default_libraries("min-surface")


add_custom_command(
        TARGET min-surface POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/nlin.ini
                ${CMAKE_CURRENT_BINARY_DIR})
