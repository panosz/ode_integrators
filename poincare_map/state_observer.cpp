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

#include "samplingCollections.hpp"
#include "armadillo_state.hpp"
#include "input_output/text_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"
#include "input_output/hdf5_io.hpp"



/* The type of container used to hold the state vector */
namespace DS
{

    using myState = armadillo_state<6>;

    struct OneForm {
        double p;
        double q;
        OneForm (double P, double Q)
            : p(P), q(Q)
        { }
    };

    class UnperturbedExtendedPendulumHamiltonian {

     private:
      double M_;

     public:
      using StateType=myState;

      explicit UnperturbedExtendedPendulumHamiltonian (double M)
          : M_{M}
      { };

      double operator() (const myState& s) const
      {
        const auto& p = s[0];
        const auto& q = s[1];
        const auto& F = s[2];
        return M_ * p * p / 2 - F * cos(q);
      }

      double dp (const myState& s) const noexcept
      {
        const auto& p = s[0];
        const auto& q = s[1];
        const auto& F = s[2];
        return M_ * p;
      }

      double dq (const myState& s) const noexcept
      {
        const auto& q = s[1];
        const auto& F = s[2];
        return F * sin(q);
      }

      double dF (const myState& s) const noexcept
      {
        const auto& q = s[1];
        return -cos(q);
      }

    };

    template<typename UnperturbedHamiltonian>
    class UnperturbedDynamicSystem {

     private:
      UnperturbedHamiltonian h_;
     public:
      using StateType = typename UnperturbedHamiltonian::StateType;

      explicit UnperturbedDynamicSystem (const UnperturbedHamiltonian& h)
          : h_{h}
      { }

      void operator() (const StateType& s, StateType& dsdt, const double /*t*/) const
      {
        const auto& p = s[0];

        dsdt[0] = -h_.dq(s);
        dsdt[1] = h_.dp(s);
        dsdt[2] = 0;
        dsdt[3] = h_.dF(s);
        dsdt[4] = p *dsdt[1];
        dsdt[5] = 1;
      }
    };

    template<typename UnperturbedHamiltonian>
    UnperturbedDynamicSystem<UnperturbedHamiltonian> makeUnperturbedDynamicSystem(const UnperturbedHamiltonian& h)
    {
      return UnperturbedDynamicSystem<UnperturbedHamiltonian>(h);
    }




}

void print_usage_string ()
{
  const auto usage_string = "usage: poincare_map input_filename integration_time perturpbation_amplitude q_harmonic chi_harmonic";
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

  inputOptions.integration_time = get_double_from_input(argv[2], "integration time");
  inputOptions.perturbation_amplitude = get_double_from_input(argv[3], "perturbation amplitude");
  inputOptions.q_harmonic = get_int_from_input(argv[4], "q_harmonic");
  inputOptions.chi_harmonic = get_int_from_input(argv[5], "chi_harmonic");

  return inputOptions;

}

class ActionIntegrationResult {
 private:
  double Action_;
  double theta_period_;
  double delta_phi_;
 public:
  ActionIntegrationResult (double Action, double theta_period, double delta_phi)
      : Action_{Action}, theta_period_{theta_period}, delta_phi_{delta_phi}
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
  ///TODO: Replace raw indices with variable tags

  const auto normalized_delta_Action = (last_point[4] - first_point[4]) * one_div_two_pi;
  const auto theta_period = last_point[5] - first_point[5];
  const auto delta_phi = last_point[3] - first_point[3];

  return ActionIntegrationResult{normalized_delta_Action, theta_period, delta_phi};
};
int main (int argc, char *argv[])
{

  const auto user_options = parse_input(argc, argv);

  const auto input_filename = user_options.input_filename;
  const auto integration_time = user_options.integration_time;
  const auto perturbation_amplitude = user_options.perturbation_amplitude;
  const auto q_harmonic = user_options.q_harmonic;
  const auto chi_harmonic = user_options.chi_harmonic;

  const auto init_states =
      get_state_from_file<DS::UnperturbedExtendedPendulumHamiltonian::StateType>(input_filename, 6);

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions(1e-12, 1e-12, 1e-5);

  auto my_sys = DS::makeUnperturbedDynamicSystem(DS::UnperturbedExtendedPendulumHamiltonian(1.0));

  std::cout << "Demonstrate integration of single init state\n";
  const auto init_state = init_states[0];
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
  std::cout << "init_F\t\tAction\tomega_theta\tomega_phi\tg_factor\n";

  for (const auto& s:init_states)
    {
      const auto action_result = action_integration(my_sys, s, user_options.integration_time, options);
      std::cout << s[2] << '\t' << action_result.Action() << '\t' << action_result.omega_theta()
                << '\t' << action_result.omega_phi() << '\t' << action_result.g_factor() << '\n';

    }

  return 0;
}
