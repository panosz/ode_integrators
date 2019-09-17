
add_library(hello_ext MODULE src/hello_ext.cpp)

target_link_libraries(hello_ext Boost::python37 Boost::numpy37 Python3::Python poincare_map )
set_target_properties(hello_ext PROPERTIES PREFIX "")

target_compile_options(hello_ext
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


configure_file(src/test_hello_ext.py ${CMAKE_CURRENT_BINARY_DIR}/test_hello_ext.py COPYONLY)
