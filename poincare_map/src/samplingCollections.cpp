//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#include "samplingCollections.hpp"

std::ostream& operator<< (std::ostream& os, const OrbitCrossOutput& orbitCrossOutput)
{

  os << "# initial point:\n";
  orbitCrossOutput.initial_point().raw_print(os);
  os << "# cross points:\n";
  orbitCrossOutput.cross_points().raw_print(os);

  return os;
}
