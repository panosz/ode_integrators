//
// Created by Panagiotis Zestanakis on 26/05/18.
//

#ifndef ODE_INTEGRATORS_PHASE_SPACE_HPP
#define ODE_INTEGRATORS_PHASE_SPACE_HPP

#include <armadillo>
#include <boost/operators.hpp>

template <int N>
class PhaseSpacePoint: boost::additive1< PhaseSpacePoint<N> ,
    boost::additive2< PhaseSpacePoint<N> , double ,
        boost::multiplicative2< PhaseSpacePoint<N> , double > > >
{
 private:
  arma::mat::fixed<N,2> state;
 public:


};

#endif //ODE_INTEGRATORS_PHASE_SPACE_HPP
