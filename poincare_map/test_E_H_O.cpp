//
// Created by Panagiotis Zestanakis on 18/09/18.
//
#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <filesystem>
#include <cmath>
#include <stdexcept>
#include <boost/tuple/tuple.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/range/combine.hpp>
#include <boost/math/special_functions/pow.hpp>
#include <boost/program_options.hpp>
namespace po = boost::program_options;


#include "samplingCollections.hpp"
#include "armadillo_state.hpp"
#include "input_output/text_file_io.hpp"
#include "system_and_poincare_surface.hpp"
#include "integration_utilities.hpp"
#include "input_output/hdf5_io.hpp"
#include "hamiltonian_dynamic_system.hpp"
#include "action_integration_result.hpp"
#include "action_integration.hpp"
#include "UnperturbedExtendedOscillatorHamiltonian.hpp"


#define TEST_E_H_O_NAME "test_E_H_0"
#define TEST_E_H_O_VERSION "1.0"

//Prints version info
void print_version_info()
{
  std::cout << TEST_E_H_O_NAME
            << ", version "
            << TEST_E_H_O_VERSION
            << '\n';
}

void print_description()
{
  print_version_info();
  std::cout<<'\n';
  std::cout << "This program is used to test the path integral calculation"
               " of the Hessian of the Hamiltonian in action space\n"
            << "The Hamiltonian used is that of the Extended Harmonic Oscillator\n"
            << "\t\t H = M p^2 + F q^2\n"
            << "Its functional form in Action space is"
            << "\t\t K = 2 Sqrt(M F) J\n";
}


// Return a group of options that will be
// allowed only on command line
po::options_description generic_options(std::string& config_file)
{
    po::options_description generic("Generic options");
    generic.add_options()
      ("version,v",
       "print version string")
      ("help,h",
       "produce help message")
      ("config,c",
        po::value<std::string>(&config_file),
        "configuration file")
      ;
    return generic;
}

std::ostream& report_flag(std::ostream& os,
                          bool flag,
                          const std::string flag_description)
{
  return os << (flag ? "Running\t": "NOT running")
            << "\t\"" << flag_description << "\"\n";
}

struct TestFlags
{
    bool specific_times=true;
    bool complete_orbit=true;
    bool complete_extended_orbit=true;
    bool analytical=true;
};

std::ostream& operator<<(std::ostream& os, const TestFlags& tf)
{
    report_flag(os,
                tf.specific_times,
                "calculation at specific times");
    report_flag(os,
                tf.complete_orbit,
                "calculate and print out a complete closed orbit");
    report_flag(os,
                tf.complete_extended_orbit,
                "calculate and print out a complete closed orbit"
                "in extended phase space");
    report_flag(os,
                tf.analytical,
                "test the numerical action integration results"
                "versus the analytical calculations");
    return os;
}

struct UserOptions {
    TestFlags test_flags{};
    std::string init_filename{};
    double integration_time{};
    double mass{};
    unsigned initial_point_index{};
};

std::ostream& operator<<(std::ostream& os, const UserOptions& uo)
{
    os<<uo.test_flags
      <<'\n'
      << "Initial points file is: "
      << uo.init_filename << '\n'
      << "Integration time is " << uo.integration_time << '\n'
      << "System mass is" << uo.mass << '\n'
      << "The initial point used in the single orbits test"
      << " is the point with index " << uo.initial_point_index << '\n';
    return os;
}

// allowed both on command line and in
// config file
po::options_description configuration_options(UserOptions& uo)

{
    po::options_description config("Configuration");
    config.add_options()
        ("specific_times",
         po::value<bool>(&uo.test_flags.specific_times)->default_value(true),
         "When on, calculate and print out an orbit at specific times")
        ("complete_orbit",
         po::value<bool>(&uo.test_flags.complete_orbit)->default_value(true),
         "When on, calculate and print out a complete closed orbit")
        ("complete_extended_orbit",
         po::value<bool>(&uo.test_flags.complete_extended_orbit)->default_value(true),
         "When on, calculate and print out a complete closed orbit"
         " in extended phase space")
        ("analytical",
         po::value<bool>(&uo.test_flags.analytical)->default_value(true),
         "When on, test the numerical action integration results"
         "versus the analytical calculations")
        ("integration_time,t",
          po::value<double>(&uo.integration_time)->required(),
          "maximum integration time")
        ("mass",
          po::value<double>(&uo.mass)->required(),
          "system mass")
        ("init,i",
          po::value<std::string>(&uo.init_filename)->required(),
          "initial positions file")
        ("index",
          po::value<unsigned>(&uo.initial_point_index)->default_value(0),
          "determines the initial point used for the single orbit tests")
        ;

    return config;
}

class TestRegister
{
public:
  explicit TestRegister(std::string name):name_{name}{};

  TestRegister& operator +=(bool test_result) noexcept
  {
    number_of_tests_++;
    if (test_result)
      passed_tests_++;
    else
      failed_tests_++;

    return *this;
  }

  unsigned total_tests() const noexcept
  {
    return number_of_tests_;
  };

  unsigned passed() const noexcept
  {
    return passed_tests_;
  };

  unsigned failed() const noexcept
  {
    return failed_tests_;
  }


  virtual ~TestRegister()
  {

      std::cout << name_<< " TEST RESULTS\n";
      std::cout << "-------------------------------------------------------------\n";
      std::cout << "*   TOTAL TESTS: " << total_tests() << '\n';
      std::cout << "*   PASSED: " << passed() << '\n';
      std::cout << "*   FAILED: " << failed() << '\n';
      std::cout << "-------------------------------------------------------------\n";
  }


private:
  std::string name_{};
  unsigned number_of_tests_ = 0;
  unsigned passed_tests_ = 0;
  unsigned failed_tests_ = 0;
};


bool double_near (double x, double y, double abs_tolerance)
{
  return std::abs(x - y) < abs_tolerance;
}

bool position_absolute_near(const DS::PhaseSpaceState& p1,
                   const DS::PhaseSpaceState& p2,
                   double tolerance)
{
  using boost::numeric::odeint::vector_space_norm_inf;
  const DS::PhaseSpaceState dp = p1-p2;
  return double_near(vector_space_norm_inf<DS::PhaseSpaceState>()(dp), 0, tolerance);
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

        UserOptions uo;
        std::string config_file;

    try {
        // Declare a group of options that will be
        // allowed only on command line
        auto generic = generic_options(config_file);

        // Declare a group of options that will be
        // allowed both on command line and in
        // config file
        auto config = configuration_options(uo);

        po::options_description all_opt_description("Allowed options");
        all_opt_description.add(generic).add(config);

        po::variables_map vm;
        po::parsed_options generic_parced=po::command_line_parser(argc, argv).
              options(generic).allow_unregistered().run();
        store(generic_parced, vm);
        notify(vm);

        if (vm.count("help"))
        {
            print_description();
            std::cout << all_opt_description << '\n';
            return 0;
        }

        if (vm.count("version"))
        {
            print_version_info();
            return 0;
        }

        std::vector<std::string> to_pass_further =
          po::collect_unrecognized(generic_parced.options, po::include_positional);

        store(po::command_line_parser(to_pass_further).
              options(config).run(), vm);

        if (vm.count("config"))
        {
          std::ifstream ifs(config_file.c_str());
            if (!ifs)
            {
              std::cerr << "Error: Can not open config file: " << config_file << '\n';
              std::cout << all_opt_description << '\n';
              return 0;
            }
            else
            {
              store(parse_config_file(ifs, config), vm);
            }
        }
        notify(vm);

        std::cout<<uo<<'\n';


    }
    catch(std::exception& e)
    {
        std::cout << e.what() << '\n';
        return 1;
    }


  // const auto user_options = parse_arguments(argc, argv);


  const auto init_states =
      get_state_from_file<DS::PhaseSpaceState>(uo.init_filename.c_str(), 4);

  if (init_states.empty())
    throw(std::runtime_error("init filename does not contain any initial states"));

  if (uo.initial_point_index >= init_states.size())
    throw(std::runtime_error("Option index greater than the number"
          " of initial states in init filename"));

  const auto init_state = init_states[uo.initial_point_index];

  std::cout << "Init States:\n" << init_states << '\n';

  const auto options = IntegrationOptions{1e-12, 1e-12, 1e-5};

  const auto myHam = DS::UnperturbedExtendedOscillatorHamiltonian(uo.mass);

  TestRegister test_register("Test Register");

  if (uo.test_flags.specific_times)
  {
    auto my_phase_space_sys = DS::makeUnperturbedDynamicSystem(myHam);
    std::cout << "\nDemonstrate integration at specific times\n";
    std::cout << "-----------------------------------------\n";


    std::cout << "init_state = " << init_state;

    const std::vector<double> times{0, 0.01,0.02,1.03,3.04};
    const auto points_at_times = orbit_points_at_times(
                                                my_phase_space_sys,
                                                init_state,
                                                times,
                                                options);

    for (const auto & point : points_at_times)
      std::cout<<point<<'\n';

    std::vector<DS::PhaseSpaceState> analytic_points_at_times{};
    for (const auto& t : times)
      analytic_points_at_times.push_back(myHam.propagate(init_state, t));


    std::cout<< "analytic calculation of positions at specific times:"<<'\n';
    for (const auto & point : analytic_points_at_times)
      std::cout<<point<<'\n';


    bool all_accurate_enough = true;
    std::cout<< "compare with analytic results:"<<'\n';
    for (const auto & zipped : boost::combine(points_at_times,
                                              analytic_points_at_times))
    {
      DS::PhaseSpaceState p_numeric;
      DS::PhaseSpaceState p_exact;
        boost::tie(p_numeric, p_exact) = zipped;

        std::cout<<"Exact position: " << p_exact;
        std::cout<<"Distance : " << p_numeric - p_exact;
        bool accurate_enough = position_absolute_near(p_numeric,p_exact,1e-8);
        std::cout << "Numeric position is"
                  <<(accurate_enough?" ":" NOT ")
                  << "accurate enough" << '\n';
        all_accurate_enough = all_accurate_enough && accurate_enough;
    }

    std::cout<< "TEST "<< (all_accurate_enough?"PASSED":"FAILED!")<<'\n';
    test_register += all_accurate_enough;




    std::cout << "\nEnd of Demonstrate integration at specific times\n";
    std::cout << "------------------------------------------------\n";
  }

  if(uo.test_flags.complete_orbit)
  {
    auto my_phase_space_sys = DS::makeUnperturbedDynamicSystem(myHam);
    std::cout << "\nDemonstrate integration in single closed orbit\n";
    std::cout << "----------------------------------------------\n";

    std::cout << "init_state = " << init_state;

    const auto closed_orbit = integrate_along_closed_orbit(
                                                my_phase_space_sys,
                                                init_state,
                                                uo.integration_time,
                                                options);

    const auto last_point = last_point_on_closed_orbit(
                                                my_phase_space_sys,
                                                init_state,
                                                uo.integration_time,
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

  if(uo.test_flags.complete_extended_orbit)
  {
    auto my_sys = DS::makeActionDynamicSystem(myHam);
    std::cout << "\nDemonstrate extended integration in single closed orbit\n";
    std::cout << "-------------------------------------------------------\n";

    const auto extended_init_state = DS::phase_to_extended_space_state(init_state);
    std::cout << "extended_init_state = " << extended_init_state;

    const auto closed_orbit = integrate_along_closed_orbit(
                                                my_sys,
                                                extended_init_state,
                                                uo.integration_time,
                                                options);

    const auto last_point = last_point_on_closed_orbit(
                                                my_sys,
                                                extended_init_state,
                                                uo.integration_time,
                                                options);

    for (const auto& s : closed_orbit)
      std::cout << s;

    std::cout << "extended_init_state " << extended_init_state
              << "extended_final_state" << closed_orbit.back() ;

    std::cout <<
      "\nEnd of Demonstrate extended integration in single closed orbit\n";
    std::cout <<
      "--------------------------------------------------------------\n";

  }

  if(uo.test_flags.analytical)
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
                           uo.integration_time,
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

    for (const auto& s:init_states)
    {
      const auto action_result = action_integration(my_sys, s, uo.integration_time, options);
      test_register += test_integration_result(s, action_result, myHam, 1e-10);

    }
 }

  return 0;
}
