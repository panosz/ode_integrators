//
// Created by Panagiotis Zestanakis on 18/09/18.
//

#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <cmath>

#include "orbit_classes.hpp"
#include "state_wrapper.hpp"
#include "state_specific_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"



/* The type of container used to hold the state vector */
namespace DS
{

    StateWrapper unperturbedExtendedHOderivatives (const StateWrapper& s)
    {
      const auto& q = s[1];
      const auto& F = s[2];

      return StateWrapper{
          -F * sin(q),
          s[0],
          0,
          -cos(q)
      };
    }

    /// \brief Models a perturbation $\Delta H = \sin(m q + n chi)$
    struct SinePerturbation {
        int q_harmonic;
        int chi_harmonic;
        SinePerturbation (int qHarmonic, int chiHarmonic)
            : q_harmonic(qHarmonic), chi_harmonic(chiHarmonic)
        { }

        double phase (const StateWrapper& s) const
        {
          const auto& q = s[1];
          const auto& chi = s[3];

          return q_harmonic * q + chi_harmonic * chi;
        }

        StateWrapper operator() (const StateWrapper& s, double /*t*/) const
        {
          const auto myCos = cos(phase(s));

          const auto dpdt = -q_harmonic * myCos;
          const auto dFdt = -chi_harmonic * myCos;

          return StateWrapper{dpdt, 0, dFdt, 0};

        }

    };

    class ExtendedHarmonicOscillator {

      double amplitude_ = 0;
      SinePerturbation perturb_;
     public:
      using StateType = StateWrapper;

      ExtendedHarmonicOscillator (double amplitude, const SinePerturbation& perturb)
          : amplitude_(amplitude), perturb_(perturb)
      { }

      void operator() (const StateType& s, StateType& dsdt, const double t) const
      {
        dsdt = unperturbedExtendedHOderivatives(s) + amplitude_ * perturb_(s, t);
      }

    };

}
//void write_to_files(const char* filename_)


void print_usage_string ()
{
  const auto usage_string = "usage: state_observer input_filename integration_time perturpbation_amplitude q_harmonic chi_harmonic";
  std::cout << usage_string << std::endl;
}

struct InputOptions {
    char *input_filename;
    double integration_time;
    double perturbation_amplitude;
    int q_harmonic;
    int chi_harmonic;

};

double get_double_from_input (const std::string& input, const std::string& parameter)
{
  double out;
  try
    {
      out = std::stod(input);
    }
  catch (...)
    {
      std::cerr << "Invalid" + parameter + ": " + input << std::endl;
      throw;
    }

  return out;
}


int get_int_from_input (const std::string& input, const std::string& parameter)
{
  int out;
  try
    {
      out = std::stoi(input);
    }
  catch (...)
    {
      std::cerr << "Invalid" + parameter + ": " + input << std::endl;
      throw;
    }

  return out;
}


InputOptions parse_input (int argc, char *argv[])
{
  InputOptions inputOptions{};

  if (argc < 6)
    print_usage_string();

  switch (argc)
    {
  case 1:
    throw std::runtime_error("input_filename must be specified.");

  case 2:
    throw std::runtime_error("integration time must be specified");

  case 3:
    throw std::runtime_error("perturbation_amplitude must be specivied");

  case 4:
    throw std::runtime_error("q_harmonic must be specivied");

  case 5:
    throw std::runtime_error("chi_harmonic must be specivied");
    }

  inputOptions.input_filename = argv[1];

  inputOptions.integration_time = get_double_from_input(argv[2],"integration time");
  inputOptions.perturbation_amplitude = get_double_from_input(argv[3],"perturbation amplitude");
  inputOptions.q_harmonic = get_int_from_input(argv[4],"q_harmonic");
  inputOptions.chi_harmonic = get_int_from_input(argv[5],"chi_harmonic");


  return inputOptions;

}

int main (int argc, char *argv[])
{

  const auto user_options = parse_input(argc, argv);

  const auto input_filename = user_options.input_filename;
  const auto integration_time = user_options.integration_time;
  const auto perturbation_amplitude = user_options.perturbation_amplitude;
  const auto q_harmonic = user_options.q_harmonic;
  const auto chi_harmonic = user_options.chi_harmonic;

  const auto init_states =
      get_state_from_file<DS::ExtendedHarmonicOscillator::StateType>(input_filename, 4);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const int zero_cross_positive_direction = 1;
  const double zero_cross_position = 0.0;

  const auto perturbation_form = DS::SinePerturbation(q_harmonic, chi_harmonic);

  const auto my_sys = DS::ExtendedHarmonicOscillator(perturbation_amplitude, perturbation_form);

  const auto my_poincare_surface = Surface{VariableTag::q,
                                           zero_cross_position,
                                           zero_cross_positive_direction};

  auto my_system_and_pc = make_system_and_poincare_surface(my_sys, my_poincare_surface);

  auto poincare_points = trace_on_poincare_surface(my_system_and_pc, init_states, integration_time, options);

  write_to_files("cross", poincare_points);

  return 0;
}
