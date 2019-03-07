cmake_minimum_required(VERSION 3.10)
project(ode_integrators)

set(CMAKE_CXX_STANDARD 17)

find_package(Boost REQUIRED)

find_package(gnuplot_iostream REQUIRED)

find_package(myUtilities REQUIRED)
find_package(Armadillo REQUIRED)


add_subdirectory(stand_alone_odeint_examples)

add_subdirectory(poincare_map)
