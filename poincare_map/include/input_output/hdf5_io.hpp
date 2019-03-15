//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#ifndef ODE_INTEGRATORS_HDF5_IO_HPP_HPP
#define ODE_INTEGRATORS_HDF5_IO_HPP_HPP

#include "prepare.hpp"
#include "armadillo_state.hpp"



std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state>& states);


void write_to_hdf5_files (const char *filename,
                          const std::vector<OrbitCrossOutput>& crossOutputs);

#endif //ODE_INTEGRATORS_HDF5_IO_HPP_HPP
