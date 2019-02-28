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
namespace DS
{
    class StateWrapper {
      using T=std::array<double,3>;
     private:
      T s_{};

     public:
      using value_type = typename T::value_type;
      using  iterator =  typename T::iterator ;
      using const_iterator = typename T::const_iterator ;

      StateWrapper (value_type x, value_type y, value_type z)
          : s_{x, y, z}
      { };
      StateWrapper () = default;

      value_type& operator[] (unsigned i)
      { return s_[i]; };

      value_type operator[] (unsigned i) const
      { return s_[i]; };

      auto begin ()
      {
        return s_.begin();
      }

      auto begin () const
      {
        return s_.begin();
      }

      auto end ()
      {
        return s_.end();
      }

      auto end () const
      {
        return s_.end();
      }

    };

    std::ostream& operator<< (std::ostream& os, const StateWrapper& state)
    {
      {
        boost::copy(state, std::ostream_iterator<typename StateWrapper::value_type>(os, " "));
        return os;
      }
    }

    class ExtendedHarmonicOscillator {
     public:
      using StateType = StateWrapper;
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

}
//void write_to_files(const char* filename_)


int main (int argc, char* argv[])
{

  if (argc<2)
    throw std::runtime_error("input_filename must be specified.");



  const auto input_filename = argv[1];

  const auto init_states =
      get_state_from_file<DS::ExtendedHarmonicOscillator::StateType >(input_filename, 3);

  std::cout <<"Init States:\n"<<init_states<<'\n';


  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const int zero_cross_positive_direction = 1;
  const double zero_cross_position = 0.0;
  const double integration_time = 10030;





  const double F_factor = 1.0;
  auto my_sys = DS::ExtendedHarmonicOscillator(F_factor);

  const auto my_poincare_surface = Surface{VariableTag::q,
                                           zero_cross_position,
                                           zero_cross_positive_direction};

  auto my_system_and_pc = make_system_and_poincare_surface(my_sys, my_poincare_surface);

  auto poincare_points = trace_on_poincare_surface(my_system_and_pc, init_states, integration_time, options);


  write_to_files("cross", poincare_points);

  return 0;
}
