//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_ORBIT_CLASSES_HPP
#define ODE_INTEGRATORS_ORBIT_CLASSES_HPP

#include <array>
#include <iostream>
#include <vector>
#include <boost/range/algorithm/copy.hpp>

namespace {
    template<typename State>
    std::ostream& operator<<(std::ostream& os, const std::vector<State> & state)
    {
      {
        boost::copy(state, std::ostream_iterator< State>(os, "\n"));
        return os;
      }
    }
}


template<typename State>
struct OrbitCrossOutput {
    State initial_point{};
    std::vector<State> cross_points{};
};



template<typename State>
std::ostream& operator<< (std::ostream& os, const OrbitCrossOutput<State>& orbitCrossOutput)
{
  os << "# initial point:\n";
  os << orbitCrossOutput.initial_point;
  os << '\n';
  os << "# cross points:\n";
  os << orbitCrossOutput.cross_points;

  return os;
}

template<typename State>
struct ApproximateAndExactCrossOutput {
    OrbitCrossOutput<State> approximate;
    OrbitCrossOutput<State> exact;
};



#endif //ODE_INTEGRATORS_ORBIT_CLASSES_HPP
