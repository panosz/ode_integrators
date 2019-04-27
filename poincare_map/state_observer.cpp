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

struct SpecialIntegrals {
    double beta = 0;
    double gamma = 0;
    double beta1 = 0;
    double beta2 = 0;
    double gamma1 = 0;
    double gamma2 = 0;
};

class ActionIntegrationResult {
 private:
  double Action_;
  double theta_period_;
  double delta_phi_;
  SpecialIntegrals integrals_;

  double dT_dJ () const
  {
    //\[Omega]/(2 \[Pi]) beta1Integral;

    return integrals_.beta1 / theta_period();
  }

  double dT_dF () const
  {
    using boost::math::double_constants::one_div_two_pi;
    //dTdF = 1/(2 \[Pi]) beta2Integral + myG*dTdJ;
    const auto g = g_factor();
    const auto beta2 = integrals_.beta2;

    return one_div_two_pi * beta2 + g * dT_dJ();
  }

 public:
  ActionIntegrationResult (double Action, double theta_period, double delta_phi, SpecialIntegrals integrals)
      : Action_{Action}, theta_period_{theta_period}, delta_phi_{delta_phi}, integrals_{integrals}
  { };
  double Action () const noexcept
  {
    return Action_;
  };
  double theta_period () const noexcept
  {
    return theta_period_;
  };
  double delta_phi () const noexcept
  {
    return delta_phi_;
  };

  double omega_theta () const
  {
    using boost::math::double_constants::two_pi;
    return two_pi / theta_period();
  };

  double g_factor () const noexcept
  {
    using boost::math::double_constants::one_div_two_pi;
    return delta_phi() * one_div_two_pi;
  };

  double omega_phi () const
  {
    return delta_phi() / theta_period();
  };

  double domega_dJ () const
  {
    using boost::math::double_constants::two_pi;
    using boost::math::pow;
    //d\[Omega]dJ = -((2 \[Pi])/T^2) dTdJ;
    return -two_pi * dT_dJ() / pow<2>(theta_period());
  }

  double d2K_dJ2 () const
  {
    using boost::math::double_constants::two_pi;
    using boost::math::pow;

    return -pow<3>(omega_theta()) * integrals_.beta1 / pow<2>(two_pi);
  }

  double one_div_two_pi_gamma() const
  {
    using boost::math::double_constants::one_div_two_pi;

    return one_div_two_pi * integrals_.gamma;
  }

  double domega_dF () const
  {
    using boost::math::double_constants::two_pi;
    using boost::math::pow;
    //d\[Omega]dF = -((2 \[Pi])/T^2) dTdF;

    return -two_pi / pow<2>(theta_period()) * dT_dF();
  }

  double d2K_dJdF () const
  {
    using boost::math::pow;
    using boost::math::double_constants::one_div_two_pi;

    const auto beta1 = integrals_.beta1;
    const auto gamma1 = integrals_.gamma1;
    const auto G = g_factor();
    const auto T = theta_period();
    const auto omega_sq = pow<2>(omega_theta());

    return omega_sq * one_div_two_pi * (gamma1 - G / T * beta1);
  }

};

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
};
int main (int argc, char *argv[])
{

  const auto user_options = parse_input(argc, argv);

  const auto input_filename = user_options.input_filename;
  const auto integration_time = user_options.integration_time;
  const auto perturbation_amplitude = user_options.perturbation_amplitude;
  const auto q_harmonic = user_options.q_harmonic;
  const auto phi_harmonic = user_options.phi_harmonic;

  const auto init_states =
      get_state_from_file<DS::UnperturbedExtendedPendulumHamiltonian::StateType>(input_filename, 12);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  auto my_sys = DS::makeUnperturbedDynamicSystem(DS::UnperturbedExtendedPendulumHamiltonian(1.0));

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
  std::cout << "init_F\t\tAction\tomega_theta\tomega_phi\tg_factor\tgamma_over_two_pi\tdomega_dJ\td2K_dJ2\tdomega_dF\td2K_dJdF\n";

  for (const auto& s:init_states)
    {
      const auto action_result = action_integration(my_sys, s, user_options.integration_time, options);
      std::cout << s[2] << '\t' << action_result.Action() << '\t' << action_result.omega_theta()
                << '\t' << action_result.omega_phi() << '\t' << action_result.g_factor()<<'\t' << action_result.one_div_two_pi_gamma()
                << '\t' << action_result.domega_dJ() << '\t' << action_result.d2K_dJ2()
                << '\t' << action_result.domega_dF()
                << '\t' << action_result.d2K_dJdF() << '\n';

    }

  return 0;
}
