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

#include <boost/operators.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include "general_file_io.hpp"
#include "orbit_classes.hpp"
#include "state_specific_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"



/* The type of container used to hold the state vector */
namespace DS
{
    class StateWrapper : boost::additive<StateWrapper, boost::additive<StateWrapper, double,
        boost::multiplicative<StateWrapper, double> > > {
      using T=std::array<double, 3>;
     private:
      T s_{0, 0, 0};

     public:
      using value_type = typename T::value_type;
      using iterator =  typename T::iterator;
      using const_iterator = typename T::const_iterator;

      StateWrapper (value_type x, value_type y, value_type z)
          : s_{x, y, z}
      { };
      StateWrapper ()
          : s_{{0, 0, 0}}
      { };

      value_type& operator[] (unsigned i)
      { return s_[i]; };

      value_type operator[] (unsigned i) const
      { return s_[i]; };

      auto begin ()
      {
        return std::begin(s_);
      }

      auto begin () const
      {
        return std::cbegin(s_);
      }

      auto end ()
      {
        return std::end(s_);
      }

      auto end () const
      {
        return std::cend(s_);
      }

      StateWrapper& operator+= (double d)
      {
        for (auto& x: s_)
          x += d;
        return *this;
      }

      StateWrapper& operator*= (double d)
      {
        for (auto& x: s_)
          x *= d;
        return *this;
      }

      StateWrapper& operator+= (const StateWrapper& other)
      {
        for (int i = 0; i < std::size(s_); ++i)
          s_[i] += other.s_[i];
        return *this;
      }

      StateWrapper operator/ (const StateWrapper& other) const
      {
        StateWrapper output{*this};
        for (unsigned long i = 0; i < std::size(s_); ++i)
          output.s_[i] /= other.s_[i];
        return output;
      }

      StateWrapper abs () const
      {
        StateWrapper output{*this};
        for (auto& x: output.s_)
          x = std::abs(x);
        return output;
      }

      double norm_inf () const
      {
        const auto s_abs = abs();
        auto max_it = boost::range::max_element(s_abs.s_);
        return *max_it;
      }

    };

    StateWrapper abs (const StateWrapper& s)
    {
      return s.abs();
    }



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
      double F_factor;
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

namespace boost { namespace numeric { namespace odeint {
            template<>
            struct vector_space_norm_inf< DS::StateWrapper >
            {
                typedef double result_type;
                double operator()( const DS::StateWrapper &s ) const
                {
                  return s.norm_inf();
                }
            };
        } } }


int main (int argc, char *argv[])
{

  if (argc < 2)
    throw std::runtime_error("input_filename must be specified.");

  const auto input_filename = argv[1];

  const auto init_states =
      get_state_from_file<DS::ExtendedHarmonicOscillator::StateType>(input_filename, 3);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const int zero_cross_positive_direction = 1;
  const double zero_cross_position = 0.0;
  const double integration_time = 10003;

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
