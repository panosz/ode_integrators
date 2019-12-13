add_custom_target(pure_python_modules)
add_dependencies(poincare_map pure_python_modules)
add_custom_command(TARGET pure_python_modules POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                   ${CMAKE_CURRENT_SOURCE_DIR}/orbital
                   ${CMAKE_CURRENT_BINARY_DIR}/orbital)
add_subdirectory(core)
