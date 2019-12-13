add_custom_target(dynamic_analysis_package)
add_dependencies(poincare_map dynamic_analysis_package)

add_custom_command(
  TARGET dynamic_analysis_package POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/dynamic_analysis/__init__.py 
                ${CMAKE_CURRENT_BINARY_DIR}/dynamic_analysis/__init__.py)


add_custom_command(
  TARGET dynamic_analysis_package POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy
                ${CMAKE_CURRENT_SOURCE_DIR}/setup.py  
                ${CMAKE_CURRENT_BINARY_DIR}/setup.py)

add_subdirectory(dynamic_analysis)
