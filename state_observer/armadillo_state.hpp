//
// Created by Panagiotis Zestanakis on 05/03/19.
//

#ifndef ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#define ODE_INTEGRATORS_ARMADILLO_STATE_HPP
#include <armadillo>
#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>


namespace DS
{
    using armadillo_state = arma::rowvec4;

}


namespace boost { namespace numeric { namespace odeint {
            template<>
            struct vector_space_norm_inf< DS::armadillo_state >
            {
                typedef double result_type;
                double operator()( const DS::armadillo_state &p ) const
                {
                  return arma::norm(p,"inf");
                }
            };
        } } }



#endif //ODE_INTEGRATORS_ARMADILLO_STATE_HPP
