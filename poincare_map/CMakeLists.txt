find_package (Python3 COMPONENTS Development)
set(Boost_ROOT $ENV{HOME}/local/boost_1_69_0)
find_package(Boost REQUIRED python37 numpy37 system program_options)

find_package(myUtilities REQUIRED)
#find armadillo
find_package(Armadillo REQUIRED)
include_directories(SYSTEM ${ARMADILLO_INCLUDE_DIRS} )

configure_file ( "${CMAKE_CURRENT_SOURCE_DIR}/cmake/version.hpp.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/include/version.hpp")

add_library(poincare_map SHARED include/version.hpp 
            src/input_output/prepare.cpp
            include/input_output/prepare.hpp 
            src/input_output/text_file_io.cpp
            include/input_output/text_file_io.hpp 
            src/samplingCollections.cpp 
            include/samplingCollections.hpp
            src/system_and_poincare_surface.cpp 
            include/system_and_poincare_surface.hpp
            src/integration_utilities.cpp 
            include/integration_utilities.hpp
            include/armadillo_state.hpp 
            src/armadillo_state.cpp 
            src/input_output/hdf5_io.cpp 
            include/input_output/hdf5_io.hpp
            src/hamiltonian_dynamic_system.cpp 
            include/hamiltonian_dynamic_system.hpp 
            src/fields_and_brackets.cpp
            include/fields_and_brackets.hpp 
            include/action_integration_result.hpp 
            src/action_integration_result.cpp
            include/action_integration.hpp 
            src/action_integration.cpp
            src/integration_options.cpp
            src/orbit_closing_check.cpp
            src/phase_space_types.cpp
            src/UnperturbedExtendedOscillatorHamiltonian.cpp
            src/UnperturbedExtendedPendulumHamiltonian.cpp)

target_compile_options(poincare_map PRIVATE
                      -Wall 
                      -Wextra 
                      -Wshadow 
                      -Wnon-virtual-dtor 
                      -pedantic
                      -Weffc++ 
                      -Wold-style-cast 
                      -Wcast-align 
                      -Wunused 
                      -Woverloaded-virtual
                      -Wconversion -Wsign-conversion 
                      -Wmisleading-indentation
                      -Wduplicated-cond 
                      -Wduplicated-branches 
                      -Wlogical-op
                      -Wnull-dereference 
                      -Wuseless-cast 
                      -Wdouble-promotion
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





add_executable(test_E_H_O test_E_H_O.cpp)
target_link_libraries(test_E_H_O PUBLIC poincare_map Boost::program_options)
target_compile_options(test_E_H_O PRIVATE
                        -Wall 
                        -Wextra 
                        -Wshadow 
                        -Wnon-virtual-dtor 
                        -pedantic
                        -Weffc++ 
                        -Wold-style-cast 
                        -Wcast-align 
                        -Wunused 
                        -Woverloaded-virtual
                        -Wconversion 
                        -Wsign-conversion 
                        -Wmisleading-indentation
                        -Wduplicated-cond 
                        -Wduplicated-branches 
                        -Wlogical-op
                        -Wnull-dereference 
                        -Wuseless-cast 
                        -Wdouble-promotion
                        -Wformat=2

                        $<$<CONFIG:Release>:
                        -O3
                        -march=native
                        >
                        )

add_subdirectory(python)
