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

    using armadillo_base_state = arma::drowvec;

    template<arma::uword N>
    using armadillo_state = armadillo_base_state::fixed<N>;

}

template<arma::uword N>
std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state<N>>& states)
{
  {

    for (const auto& s : states)
      s.raw_print(os);
    return os;
  }
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
