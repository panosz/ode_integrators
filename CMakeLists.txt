cmake_minimum_required(VERSION 3.10)
project(ode_integrators)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

find_package(gnuplot_iostream REQUIRED)

find_package(myUtilities REQUIRED)
find_package(Armadillo REQUIRED)


add_subdirectory(stand_alone_odeint_examples)


add_executable(poincare_map poincare_map/state_observer.cpp poincare_map/input_output/general_file_io.cpp poincare_map/input_output/general_file_io.hpp poincare_map/input_output/state_specific_file_io.cpp poincare_map/input_output/state_specific_file_io.hpp poincare_map/orbit_classes.cpp poincare_map/orbit_classes.hpp poincare_map/system_and_poincare_surface.cpp poincare_map/system_and_poincare_surface.hpp poincare_map/integration_utilities.cpp poincare_map/integration_utilities.hpp poincare_map/state_wrapper.cpp poincare_map/state_wrapper.hpp poincare_map/armadillo_state.hpp poincare_map/armadillo_state.cpp poincare_map/input_output/hdf5_io.cpp poincare_map/input_output/hdf5_io.hpp)
target_link_libraries(poincare_map Boost::boost myUtilities::myUtilities stdc++fs armadillo)
target_compile_options(poincare_map PRIVATE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
        -Weffc++ -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wmisleading-indentation
        -Wduplicated-cond -Wduplicated-branches -Wlogical-op
        -Wnull-dereference -Wuseless-cast -Wdouble-promotion
        -Wformat=2

        $<$<CONFIG:Release>:
        -O3
        -march=native
        >
        )

