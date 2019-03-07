//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#ifndef ODE_INTEGRATORS_HDF5_IO_HPP_HPP
#define ODE_INTEGRATORS_HDF5_IO_HPP_HPP

#include "general_file_io.hpp"
#include "armadillo_state.hpp"

void write_to_hdf5_files (const char *filename,
                          const std::vector<OrbitCrossOutput<DS::armadillo_state>>& crossOutputs);

#endif //ODE_INTEGRATORS_HDF5_IO_HPP_HPP
