add_executable("MF-stokes" MF-stokes.cc bc_extension.hh driver.hh)
target_link_dune_default_libraries("MF-stokes")


add_custom_command(
        TARGET MF-stokes POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_SOURCE_DIR}/src/prepreka1.msh
                ${CMAKE_CURRENT_BINARY_DIR})
