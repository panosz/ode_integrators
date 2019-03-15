//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#ifndef ODE_INTEGRATORS_HDF5_IO_HPP_HPP
#define ODE_INTEGRATORS_HDF5_IO_HPP_HPP

#include "prepare.hpp"
#include "armadillo_state.hpp"


arma::mat matrix_from_collection_of_armadillo_states (const std::vector<DS::armadillo_state >& states);


struct ArmaOrbitCrossOutput
{
    DS::armadillo_state initial_point{};
    arma::mat cross_points{};

    ArmaOrbitCrossOutput() = default;
    explicit ArmaOrbitCrossOutput(const OrbitCrossOutput<DS::armadillo_state > & orbitCrossOutput);

};


std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state>& states);


void write_to_hdf5_files (const char *filename,
                          const std::vector<OrbitCrossOutput<DS::armadillo_state>>& crossOutputs);

#endif //ODE_INTEGRATORS_HDF5_IO_HPP_HPP
