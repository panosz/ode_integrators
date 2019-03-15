//
// Created by Panagiotis Zestanakis on 05/03/19.
//

#ifndef ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#define ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#include <armadillo>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>
#include "samplingCollections.hpp"

namespace DS
{

    const unsigned STATE_DIMENSIONS=4;

    using armadillo_base_state = arma::drowvec;
    using armadillo_state = armadillo_base_state::fixed<STATE_DIMENSIONS>;

}

namespace boost
{
    namespace numeric
    {
        namespace odeint
        {
            template<arma::uword N>
            struct vector_space_norm_inf<arma::drowvec::fixed<N>> {
                typedef double result_type;
                double operator() (const arma::drowvec::fixed<N>& p) const
                {
                  return arma::norm(p, "inf");
                }
            };
        }
    }
}



#endif //ODE_INTEGRATORS_ARMADILLO_STATE_HPP
