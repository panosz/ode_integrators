//
// Created by Panagiotis Zestanakis on 18/09/18.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include <vector>
#include <stdexcept>
#include <utility>

#include <cmath>

#include "general_file_io.hpp"
#include "orbit_classes.hpp"
#include "state_specific_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"



/* The type of container used to hold the state vector */

class ExtendedHarmonicOscillator {
 public:
  using StateType = std::array<double,3>;
  const double F_factor;
  explicit ExtendedHarmonicOscillator (double f)
      : F_factor{f}
  { };

  void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
  {
    dsdt[0] = -F_factor * sin(s[1]);
    dsdt[1] = s[0];
    dsdt[2] = -cos(s[1]);

  }
};




//void write_to_files(const char* filename_)

int main ()
{

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const int zero_cross_positive_direction = 1;
  const double zero_cross_position = 0.0;
  const double integration_time = 1000;

  std::vector<ExtendedHarmonicOscillator::StateType >
      init_states{
      {1.,         0, 0},
      {1.09388829, 0, 0}, // G=3/4
      {1.4,        0, 0},
      {.08,        0, 0},
      {.7,         0, 0},
      {1.23941324,0,0}, // G=2/3
      {1.15,0,0},
      {0.5,0,0},
      {0.2,0,0},
      {0.1,0,0},
      {1.45956452,0,0}, //G = 1/2
      {1.5,0,0}
  };
  const double F_factor = 1.0;
  auto my_sys = ExtendedHarmonicOscillator(F_factor);

  const auto my_poincare_surface = Surface{VariableTag::q,
                                           zero_cross_position,
                                           zero_cross_positive_direction};

  auto my_system_and_pc = make_system_and_poincare_surface(my_sys, my_poincare_surface);

  auto poincare_points = trace_on_poincare_surface(my_system_and_pc, init_states, integration_time, options);

//  std::cout << "we have as many as " << poincare_points.approximate.cross_points.size() << " zero cross points\n";
//
//  std::cout << "approximate cross output:\n";
//  std::cout << poincare_points.approximate << '\n';
//
//  std::cout << "exact cross output:\n";
//  std::cout << poincare_points.exact << '\n';

  write_to_files("cross", poincare_points);

}
