//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_SAMPLING_COLLECTIONS_HPP
#define ODE_INTEGRATORS_SAMPLING_COLLECTIONS_HPP

#include <array>
#include <iostream>
#include <vector>
#include <boost/range/algorithm/copy.hpp>
#include <armadillo>

namespace
{

}

template<typename State>
arma::drowvec rowvec_from_state (const State& state)
{
  auto out = arma::drowvec(state.size());
  boost::copy(state, out.begin());
  return out;
}

template<typename State>
arma::mat matrix_from_collection_of_states (const std::vector<State>& states)
{
  arma::mat out;
  const auto number_of_states = states.size();

  if (number_of_states == 0)
    return out;

  const auto number_of_dimensions = states[0].size();

  out.set_size(number_of_states, number_of_dimensions);

  for (long long unsigned i = 0; i < number_of_states; ++i)
    {
      out.row(i) = rowvec_from_state(states[i]);
    }
  return out;
}

class OrbitCrossOutput {
 private:
  arma::drowvec initial_point_{};
  arma::mat cross_points_{};
 public:
  OrbitCrossOutput (const arma::drowvec& initial_point, const arma::mat& cross_points)
      : initial_point_(initial_point), cross_points_(cross_points)
  { }

  const arma::drowvec& initial_point () const noexcept
  {
    return initial_point_;
  }
  const arma::mat& cross_points () const noexcept
  {
    return cross_points_;
  }

};

template<typename State>
OrbitCrossOutput make_OrbitCrossOutput (const State& init_point, const std::vector<State>& cross_points)
{
  return OrbitCrossOutput(rowvec_from_state(init_point), matrix_from_collection_of_states(cross_points));
}

std::ostream& operator<< (std::ostream& os, const OrbitCrossOutput& orbitCrossOutput);

#endif //ODE_INTEGRATORS_SAMPLING_COLLECTIONS_HPP
