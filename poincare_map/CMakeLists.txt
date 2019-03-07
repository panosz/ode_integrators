
find_package(Boost REQUIRED)

find_package(myUtilities REQUIRED)
find_package(Armadillo REQUIRED)




add_executable(poincare_map state_observer.cpp input_output/general_file_io.cpp
        input_output/general_file_io.hpp input_output/state_specific_file_io.cpp
        input_output/state_specific_file_io.hpp orbit_classes.cpp orbit_classes.hpp
        system_and_poincare_surface.cpp system_and_poincare_surface.hpp
        integration_utilities.cpp integration_utilities.hpp state_wrapper.cpp state_wrapper.hpp
        armadillo_state.hpp armadillo_state.cpp input_output/hdf5_io.cpp input_output/hdf5_io.hpp)
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

