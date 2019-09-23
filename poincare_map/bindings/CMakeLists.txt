
add_library(action_integration_ext MODULE src/action_integration_ext.cpp include/ActionIntegrationResultDecorator.hpp)

target_link_libraries(action_integration_ext Boost::python37 Boost::numpy37 Python3::Python poincare_map )
set_target_properties(action_integration_ext PROPERTIES PREFIX "")

target_compile_options(action_integration_ext
        PRIVATE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
        -Weffc++ -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wmisleading-indentation
        -Wduplicated-cond
        #-Wduplicated-branches
        -Wlogical-op
        -Wnull-dereference -Wuseless-cast -Wdouble-promotion
        -Wformat=2

        $<$<CONFIG:Release>:
        -O2
        -march=native
        >
        )
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}") 
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}") 
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}") 


configure_file(src/test_action_integration_ext.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_ext.py COPYONLY)
configure_file(src/test_action_integration_numerics.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_numerics.py COPYONLY)
