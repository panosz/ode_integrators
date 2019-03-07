//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#include "armadillo_state.hpp"
arma::mat matrix_from_collection_of_armadillo_states (const std::vector<DS::armadillo_state>& states)
{
  arma::mat out;
  const auto number_of_states = states.size();
  out.set_size(number_of_states, DS::STATE_DIMENSIONS);

  for (long long unsigned i = 0; i < number_of_states; ++i)
    {
      out.row(i) = states[i];
    }
  return out;
}

std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state>& states)
{
  {

    for (const auto& s : states)
      s.raw_print(os);
    return os;
  }
}

ArmaOrbitCrossOutput::ArmaOrbitCrossOutput (const OrbitCrossOutput<DS::armadillo_state>& orbitCrossOutput)
    :
    initial_point{orbitCrossOutput.initial_point},
    cross_points{matrix_from_collection_of_armadillo_states(orbitCrossOutput.cross_points)}
{}
