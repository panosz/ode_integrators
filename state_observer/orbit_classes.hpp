//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_ORBIT_CLASSES_HPP
#define ODE_INTEGRATORS_ORBIT_CLASSES_HPP

#include <array>
#include <iostream>
#include <vector>
#include <boost/range/algorithm/copy.hpp>

template<typename State>
struct OrbitCrossOutput {
    State initial_point{};
    std::vector<State> cross_points{};
};

template<typename State>
std::ostream& output_State_to_stream(std::ostream& os, const State & state)
{
  {
    boost::copy(state, std::ostream_iterator<typename State::value_type>(os, " "));
    return os;
  }
}

template<typename State>
std::ostream& operator<< (std::ostream& os, const OrbitCrossOutput<State>& orbitCrossOutput)
{
  os << "# initial point:\n";
  output_State_to_stream(os,orbitCrossOutput.initial_point);
  os << '\n';
  os << "# cross points:\n";
  for (const auto& item : orbitCrossOutput.cross_points)
    {
      output_State_to_stream(os, item);
      os << '\n';
    }
  return os;
}

template<typename State>
struct ApproximateAndExactCrossOutput {
    OrbitCrossOutput<State> approximate;
    OrbitCrossOutput<State> exact;
};



#endif //ODE_INTEGRATORS_ORBIT_CLASSES_HPP
