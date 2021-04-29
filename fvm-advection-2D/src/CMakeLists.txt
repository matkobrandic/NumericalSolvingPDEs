add_executable("fvm" fvm.cc data.hh evolve.hh initialize.hh)
target_link_dune_default_libraries("fvm")

# Ovo kopira samo jednom
#file(COPY fvm.input rectangle.dgf DESTINATION  "${CMAKE_BINARY_DIR}/src/")

# Ne radi
add_custom_target(copy ALL COMMAND ${CMAKE_COMMAND} -E copy "${CMAKE_SOURCE_DIR}/src/fvm.input" "${CMAKE_SOURCE_DIR}/src/rectangle.dgf"  "${CMAKE_BINARY_DIR}/src/")
