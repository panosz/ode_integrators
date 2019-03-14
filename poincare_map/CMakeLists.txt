
find_package(Boost REQUIRED)

find_package(myUtilities REQUIRED)
find_package(Armadillo REQUIRED)


configure_file (
        "${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/version.hpp"
)

add_library(poincare_map SHARED include/version.hpp src/input_output/general_file_io.cpp
        include/input_output/general_file_io.hpp src/input_output/state_specific_file_io.cpp
        include/input_output/state_specific_file_io.hpp src/samplingCollections.cpp include/samplingCollections.hpp
        src/system_and_poincare_surface.cpp include/system_and_poincare_surface.hpp
        src/integration_utilities.cpp include/integration_utilities.hpp src/state_wrapper.cpp include/state_wrapper.hpp
        include/armadillo_state.hpp src/armadillo_state.cpp src/input_output/hdf5_io.cpp include/input_output/hdf5_io.hpp)

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


target_include_directories(
        poincare_map PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include/poincare_map>
        PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
)



target_link_libraries(poincare_map Boost::boost myUtilities::myUtilities stdc++fs armadillo)


install (
        TARGETS poincare_map
        EXPORT poincare_mapTargets
        ARCHIVE  DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY  DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/poincare_map)





add_executable(generalized_pendulum state_observer.cpp)
target_link_libraries(generalized_pendulum PUBLIC poincare_map)
target_compile_options(generalized_pendulum PRIVATE
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