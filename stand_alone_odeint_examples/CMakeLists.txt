find_package(Boost REQUIRED)

find_package(gnuplot_iostream REQUIRED)

find_package(myUtilities REQUIRED)
find_package(Armadillo REQUIRED)



add_executable(harmonic_oscilator harmonic_oscilator.cpp utilities.hpp)
target_link_libraries(harmonic_oscilator Boost::boost gnuplot_iostream::gnuplot_iostream)

target_compile_options(harmonic_oscilator PRIVATE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
        -Weffc++ -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wmisleading-indentation
        -Wduplicated-cond -Wduplicated-branches -Wlogical-op
        -Wnull-dereference -Wuseless-cast -Wdouble-promotion
        -Wformat=2

        $<$<CONFIG:Release>:
        -O2
        -march=native
        >
        )


add_executable(karney_system karney_system.cpp utilities.hpp)
target_link_libraries(karney_system Boost::boost gnuplot_iostream::gnuplot_iostream)

target_compile_options(karney_system PRIVATE
        -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic
        -Weffc++ -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual
        -Wconversion -Wsign-conversion -Wmisleading-indentation
        -Wduplicated-cond -Wduplicated-branches -Wlogical-op
        -Wnull-dereference -Wuseless-cast -Wdouble-promotion
        -Wformat=2

        $<$<CONFIG:Release>:
        -O2
        -march=native
        >
        )



include_directories(${ARMADILLO_INCLUDE_DIRS})
add_executable(henon_heiles henon_heilles/henon_heiles.cpp utilities.hpp
        henon_heilles/henon_heiles_impl.cpp henon_heilles/henon_heiles_impl.hpp henon_heilles/phase_space.hpp)
target_link_libraries(henon_heiles Boost::boost gnuplot_iostream::gnuplot_iostream myUtilities::myUtilities ${ARMADILLO_LIBRARIES})

target_compile_options(henon_heiles PRIVATE
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
