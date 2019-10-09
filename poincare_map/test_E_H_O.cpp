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
#include "action_integration.hpp"

void print_usage_string (const std::string& program_name)
{
  const auto usage_string =
      "usage: " + program_name + " input_filename integration_time";
  std::cout << usage_string << std::endl;
}

struct InputOptions {
    char *input_filename;
    double integration_time;
//    double perturbation_amplitude;
//    int q_harmonic;
//    int phi_harmonic;

};

double get_double_from_argument (const std::string& input, const std::string& parameter)
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

int get_int_from_argument (const std::string& input, const std::string& parameter)
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

InputOptions parse_arguments (int argc, char **argv)
{
  InputOptions inputOptions{};

  if (argc < 3)
    print_usage_string(argv[0]);

  switch (argc)
    {
  case 1:
    throw std::runtime_error("input_filename must be specified.");

  case 2:
    throw std::runtime_error("integration time must be specified");
    }

  inputOptions.input_filename = argv[1];

  inputOptions.integration_time = get_double_from_argument(argv[2], "integration time");

  return inputOptions;

}

bool double_near (double x, double y, double abs_tolerance)
{
  return std::abs(x - y) < abs_tolerance;
}

bool test_numerical_vs_analytical_result (const std::string& quantity_name,
                                      double numerical_value,
                                      double analytical_value, double tolerance)
{
  std::cout << quantity_name + ":\n";
  std::cout << "numerical_value: " << numerical_value << '\n';
  std::cout << "analytical_value: " << analytical_value << '\n';
  std::cout << "error: " << numerical_value - analytical_value << '\n';

  auto passed = double_near(analytical_value, numerical_value, tolerance);

  std::cout << (passed ? "PASSED" : "FAILED") << '\n';
  std::cout << "========================================\n";
  return passed;

}

bool
test_integration_result (const DS::PhaseSpaceState& init_state,
                         const ActionIntegrationResult& action_integration_result,
                         const DS::UnperturbedExtendedOscillatorHamiltonian& ham,
                         double tolerance)
{
  std::cout << "Begin Testing Orbit starting from :\n"
            << init_state << '\n';
  std::cout << "========================================\n";

  bool passed = true;

  passed = passed
           && test_numerical_vs_analytical_result("Action", action_integration_result.Action(), ham.action(init_state), tolerance);
  passed = passed
           && test_numerical_vs_analytical_result("omega", action_integration_result.omega_theta(), ham.dKdJ(init_state), tolerance);
  passed = passed
           && test_numerical_vs_analytical_result("omega_phi", action_integration_result.omega_phi(), ham.dKdF(init_state), tolerance);
  passed = passed
           && test_numerical_vs_analytical_result("d2KdJ2", action_integration_result.d2K_dJ2(), ham.d2KdJ2(init_state), tolerance);
  passed = passed
           && test_numerical_vs_analytical_result("d2KdF2", action_integration_result.d2K_dF2(), ham.d2KdF2(init_state), tolerance);
  passed = passed
           && test_numerical_vs_analytical_result("d2KdJdF", action_integration_result.d2K_dJdF(), ham.d2KdJdF(init_state), tolerance);

  std::cout << "Test for Orbit starting from : " << init_state << '\n';
  std::cout << (passed ? "PASSED" : "FAILED") << '\n';
  std::cout << "================================================================================\n";

  return passed;

}

int main (int argc, char *argv[])
{

  std::cout << " TEST EXTENDED HARMONIC OSCILLATOR\n"
            << "----------------------------------------\n"
            << "This program is used to test the path integral calculation of the Hessian of the Hamiltonian in action space\n"
            << "The Hamiltonian used is that of the Extended Harmonic Oscillator\n"
            << "\t\t H = M p^2 + F q^2\n"
            << "Its functional form in Action space is"
            << "\t\t K = 2 Sqrt(M F) J\n";

  const auto user_options = parse_arguments(argc, argv);

  const auto input_filename = user_options.input_filename;

  const auto init_states =
      get_state_from_file<DS::PhaseSpaceState>(input_filename, 4);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  const auto myHam = DS::UnperturbedExtendedOscillatorHamiltonian(1.4);


  { // Uncomment this block, to calculate and print out
    // a complete closed orbit
    auto my_phase_space_sys = DS::makeUnperturbedDynamicSystem(myHam);
    std::cout << "\nDemonstrate integration in single closed orbit\n";
    std::cout << "----------------------------------------------\n";

    const auto init_state = init_states[10];

    std::cout << "init_state = " << init_state;

    const auto closed_orbit = integrate_along_closed_orbit(
                                                my_phase_space_sys,
                                                init_state,
                                                user_options.integration_time,
                                                options);

    const auto last_point = last_point_on_closed_orbit(
                                                my_phase_space_sys,
                                                init_state,
                                                user_options.integration_time,
                                                options);

    for (const auto& s : closed_orbit)
      std::cout << s;

    std::cout << "init_state " << init_state
              << "final_state" << closed_orbit.back() ;

    std::cout <<
      "\nEnd of Demonstrate integration in single closed orbit\n";
    std::cout <<
      "-----------------------------------------------------\n";
  }

  { // Uncomment this block, to calculate and print out
    // a complete closed orbit in extended phase space
    auto my_sys = DS::makeActionDynamicSystem(myHam);
    std::cout << "\nDemonstrate extended integration in single closed orbit\n";
    std::cout << "-------------------------------------------------------\n";

    const auto init_state = DS::phase_to_extended_space_state(init_states[10]);
    std::cout << "init_state = " << init_state;

    const auto closed_orbit = integrate_along_closed_orbit(
                                                my_sys,
                                                init_state,
                                                user_options.integration_time,
                                                options);

    const auto last_point = last_point_on_closed_orbit(
                                                my_sys,
                                                init_state,
                                                user_options.integration_time,
                                                options);

    for (const auto& s : closed_orbit)
      std::cout << s;

    std::cout << "init_state " << init_state
              << "final_state" << closed_orbit.back() ;

    std::cout <<
      "\nEnd of Demonstrate extended integration in single closed orbit\n";
    std::cout <<
      "--------------------------------------------------------------\n";

  }

  { // Uncomment this block, to test the numerical action integration results
    // versus the analytical calculations
    std::cout << "Compare analytical and numerical calculations\n";

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


    auto my_sys = DS::makeActionDynamicSystem(myHam);

    for (const auto& s:init_states)
    {
      const auto d2KdJ2_analytical = myHam.d2KdJ2(s);
      const auto d2KdJdF_analytical = myHam.d2KdJdF(s);
      const auto omega_phi_analytical = myHam.dKdF(s);
      const auto d2KdF2_analytical = myHam.d2KdF2(s);

      const auto action_result =
        action_integration(my_sys,
                           s,
                           user_options.integration_time,
                           options);

      std::cout << s[2] << "\t\t"
        << action_result.Action() << "\t\t"
        << action_result.omega_theta() << "\t\t"
        << action_result.omega_phi() << "\t\t"
        << omega_phi_analytical << "\t\t"
        << action_result.g_factor() << "\t\t"
        << action_result.one_div_two_pi_gamma() << "\t\t"
        << action_result.d2K_dJ2() << "\t\t"
        << d2KdJ2_analytical << "\t\t"
        << action_result.domega_dF() << "\t\t"
        << action_result.d2K_dJdF() << "\t\t"
        << d2KdJdF_analytical << "\t\t"
        << action_result.d2K_dF2() << "\t\t"
        << d2KdF2_analytical << '\n';

    }

    std::cout << "\n\n\n";
    std::cout << "**************************************\n";
    std::cout << "**********     RUN TESTS    **********\n";
    std::cout << "**************************************\n";

    int no_of_passed_tests = 0;
    int no_of_failed_tests = 0;
    int no_of_tests = 0;
    for (const auto& s:init_states)
    {
      const auto action_result = action_integration(my_sys, s, user_options.integration_time, options);
      const auto passed = test_integration_result(s, action_result, myHam, 1e-10);
      no_of_passed_tests += passed;
      no_of_failed_tests += ! passed;
      no_of_tests +=1;

    }
    std::cout << "******************************************************************************************\n";
    std::cout << "*   TOTAL TESTS: " << no_of_tests << '\n';
    std::cout << "*   PASSED: " << no_of_passed_tests<< '\n';
    std::cout << "*   FAILED: " << no_of_failed_tests<< '\n';
    std::cout << "******************************************************************************************\n";
  }

  return 0;
}
