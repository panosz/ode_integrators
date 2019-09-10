//
// Created by Panagiotis Zestanakis on 18/09/18.
//
#include <vector>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <cmath>
#include <stdexcept>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/pow.hpp>

#include "samplingCollections.hpp"
#include "armadillo_state.hpp"
#include "input_output/text_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"
#include "input_output/hdf5_io.hpp"
#include "hamiltonian_dynamic_system.hpp"
#include "action_integration_result.hpp"

void print_usage_string ()
{
  const auto usage_string = "usage: poincare_map input_filename integration_time perturpbation_amplitude q_harmonic phi_harmonic";
  std::cout << usage_string << std::endl;
}

struct InputOptions {
    char *input_filename;
    double integration_time;
    double perturbation_amplitude;
    int q_harmonic;
    int phi_harmonic;

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
    throw std::runtime_error("phi_harmonic must be specivied");
    }

  inputOptions.input_filename = argv[1];

  inputOptions.integration_time = get_double_from_input(argv[2], "integration time");
  inputOptions.perturbation_amplitude = get_double_from_input(argv[3], "perturbation amplitude");
  inputOptions.q_harmonic = get_int_from_input(argv[4], "q_harmonic");
  inputOptions.phi_harmonic = get_int_from_input(argv[5], "phi_harmonic");

  return inputOptions;

}



template<typename System>
ActionIntegrationResult
action_integration (System sys,
                    const typename System::StateType& first_point,
                    double max_time,
                    IntegrationOptions integrationOptions)
{
  using boost::math::double_constants::one_div_two_pi;

  const auto last_point = last_point_on_closed_orbit(sys, first_point, max_time, integrationOptions);
  const auto delta_point = last_point - first_point;

  SpecialIntegrals integrals{};

  integrals.beta = delta_point[static_cast<unsigned>(CoordinateTag::beta)];
  integrals.beta1 = delta_point[static_cast<unsigned>(CoordinateTag::beta1)];
  integrals.beta2 = delta_point[static_cast<unsigned>(CoordinateTag::beta2)];

  integrals.gamma = delta_point[static_cast<unsigned>(CoordinateTag::gamma)];
  integrals.gamma1 = delta_point[static_cast<unsigned>(CoordinateTag::gamma1)];
  integrals.gamma2 = delta_point[static_cast<unsigned>(CoordinateTag::gamma2)];

  const auto delta_J = delta_point[static_cast<unsigned>(CoordinateTag::J)];
  const auto theta_period = delta_point[static_cast<unsigned>(CoordinateTag::t)];
  const auto delta_phi = delta_point[static_cast<unsigned>(CoordinateTag::phi)];

  const auto normalized_delta_Action = delta_J * one_div_two_pi;

  return ActionIntegrationResult{normalized_delta_Action, theta_period, delta_phi, integrals};
}

int main (int argc, char *argv[])
{

  const auto user_options = parse_input(argc, argv);

  const auto input_filename = user_options.input_filename;

  const auto init_states =
      get_state_from_file<DS::UnperturbedExtendedPendulumHamiltonian::StateType>(input_filename, 12);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const auto myHam = DS::UnperturbedExtendedOscillatorHamiltonian(1.0);
  auto my_sys = DS::makeUnperturbedDynamicSystem(myHam);

  std::cout << "Demonstrate integration of single init state\n";
  const auto init_state = init_states[10];
  std::cout << "init_state = " << init_state;

  const auto closed_orbit = integrate_along_closed_orbit(my_sys, init_state, user_options.integration_time, options);
  const auto last_point = last_point_on_closed_orbit(my_sys, init_state, user_options.integration_time, options);

  for (const auto& s : closed_orbit)
    std::cout << s;

  std::cout << "init_state " << init_state <<
            "final_state" << closed_orbit.back() <<
            "only final state" << last_point;

  std::cout << "\nEnd Demonstratie integration of single init state\n";

  std::cout << "Action integration for all states\n\n";
  std::cout << "init_F" << "\t\t"
            << "Action" << "\t\t"
            << "omega_theta" << "\t\t"
            << "omega_phi" << "\t\t"
            << "omega_phi_analytic" << "\t\t"
            << "g_factor" << "\t\t"
            << "gamma_over_two_pi" << "\t\t"
            << "d2K_dJ2" << "\t\t"
            << "domega_dJ_analytic" << "\t\t"
            << "domega_dF" << "\t\t"
            << "d2K_dJdF" << "\t\t"
            << " d2K_dJdF_analytic" << "\t\t"
            << "d2K_dF2" << "\t\t"
            << " d2K_dF2_analytic" << "\n";

  for (const auto& s:init_states)
    {
      const auto d2KdJ2_analytic = myHam.d2KdJ2(s);
      const auto d2KdJdF_analytic = myHam.d2KdJdF(s);
      const auto omega_phi_analytic = myHam.dKdF(s);
      const auto d2KdF2_analytic = myHam.d2KdF2(s);

      const auto action_result = action_integration(my_sys, s, user_options.integration_time, options);
      std::cout << s[2] << "\t\t"
                << action_result.Action() << "\t\t"
                << action_result.omega_theta() << "\t\t"
                << action_result.omega_phi() << "\t\t"
                << omega_phi_analytic << "\t\t"
                << action_result.g_factor() << "\t\t"
                << action_result.one_div_two_pi_gamma() << "\t\t"
                << action_result.d2K_dJ2() << "\t\t"
                << d2KdJ2_analytic << "\t\t"
                << action_result.domega_dF() << "\t\t"
                << action_result.d2K_dJdF() << "\t\t"
                << d2KdJdF_analytic << "\t\t"
                << action_result.d2K_dF2() << "\t\t"
                << d2KdF2_analytic << '\n';

    }

  return 0;
}
