
add_library(action_integration_ext MODULE 
            src/action_integration_ext.cpp 
            include/ActionIntegrationResultDecorator.hpp
            src/state_bindings.cpp
            include/state_bindings.hpp
            src/arma_state_bindings.cpp
            src/integration_options_bindings.cpp
            include/integration_options_bindings.hpp
            src/ActionIntegrationResultDecorator.cpp
            src/dynamic_system_bindings.cpp
            src/hamiltonians_bindings.cpp
            )

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
target_include_directories(
        action_integration_ext PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/poincare_map/bindings/include>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)
message(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}") 
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}") 
message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}") 


# configure_file(src/test_action_integration_ext.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_ext.py COPYONLY)
# configure_file(src/test_action_integration_numerics.py ${CMAKE_CURRENT_BINARY_DIR}/test_action_integration_numerics.py COPYONLY)
configure_file(__init__.py ${CMAKE_CURRENT_BINARY_DIR}/__init__.py COPYONLY)
configure_file(setup.py ${CMAKE_CURRENT_BINARY_DIR}/setup.py COPYONLY)
