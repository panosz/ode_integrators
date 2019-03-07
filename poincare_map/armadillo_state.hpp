//
// Created by Panagiotis Zestanakis on 05/03/19.
//

#ifndef ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#define ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#include <armadillo>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include "orbit_classes.hpp"

namespace DS
{

    const unsigned STATE_DIMENSIONS=4;

    using armadillo_state = arma::drowvec::fixed<STATE_DIMENSIONS>;

}

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            template<>
            struct vector_space_norm_inf<DS::armadillo_state> {
                typedef double result_type;
                double operator() (const DS::armadillo_state& p) const
                {
                  return arma::norm(p, "inf");
                }
            };
        }
    }
}


arma::mat matrix_from_collection_of_armadillo_states (const std::vector<DS::armadillo_state >& states);


struct ArmaOrbitCrossOutput
{
    DS::armadillo_state initial_point{};
    arma::mat cross_points{};

    ArmaOrbitCrossOutput() = default;
    explicit ArmaOrbitCrossOutput(const OrbitCrossOutput<DS::armadillo_state > & orbitCrossOutput);

};


std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state>& states);

#endif //ODE_INTEGRATORS_ARMADILLO_STATE_HPP
