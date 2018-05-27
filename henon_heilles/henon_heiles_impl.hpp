//
// Created by Panagiotis Zestanakis on 25/05/18.
//

#ifndef ODE_INTEGRATORS_HENON_HEILES_IMPL_HPP
#define ODE_INTEGRATORS_HENON_HEILES_IMPL_HPP

#include <iostream>
#include <array>
#include <vector>
#include "../utilities.hpp"

using vector_type = std::array<double, 2>;


class Henon_Heiles {

  static double h0 (const vector_type& q, const vector_type& p);

  static double h1 (const vector_type& q);

  static vector_type dh0dq (const vector_type& q);

  static vector_type dh0dp (const vector_type& p);

  static vector_type dh1dq (const vector_type& q);

 public :

  static double Hamiltonian (const vector_type& q, const vector_type& p);

  static vector_type dHdp (const vector_type& q, const vector_type& p);

  static void dpdt (const vector_type& q, vector_type& p_time_derivative);

  static void dqdt (const vector_type& p, vector_type& q_time_derivative);

};

double kinetic_energy_at_y (double total_energy, double y);

double derivative_of_kinetic_energy_at_y (double y);

std::pair<double, double> kinetic_energy_and_derivative (double total_energy, double y);

double kinetic_energy_root (double total_energy, double yguess);

using State = std::array<double, 4>;

vector_type get_q(const State & state);


vector_type get_p(const State & state);


std::vector<State>
make_initial_y_px_at_E (int numberOfPositions, double energy_level);;

std::ostream& operator<< (std::ostream& os, const vector_type& v);

std::ostream& operator<< (std::ostream& os, const State& s);

class TrajectoryDerivatives : public Henon_Heiles {
 public:
  virtual void operator() (const State& s, State& dsdt, const double /* t */);

};

class CrossingDerivatives : public TrajectoryDerivatives
{
 public:
  void operator() (const State& s, State& dsdt, const double t) override;

};




class SurfaceCrossEventObserver : public PanosOde::Utils::PoincareSurfaceCrossEventObserver<State> {
 protected:
  double distance_from_surface (const State& state) const override;

  bool event_predicate (double distance) const override;
};

std::pair<std::vector<double>,std::vector<double> >
henon_heiles_poincare_surface( double total_energy, double time_integration, size_t noOfInitialPoints);


#endif //ODE_INTEGRATORS_HENON_HEILES_IMPL_HPP
