//
// Created by Panagiotis Zestanakis on 24/05/18.
//
#include "henon_heiles_impl.hpp"
#include <iostream>
#include "gnuplot-iostream.h"
#include <array>

void
plot_henon_heiles_poincare_surface (Gnuplot& g1, float energy_denominator, int noOfInitialPoints, double time_integration)
{
  const double total_energy = 1. / energy_denominator;
  auto[y, py] = henon_heiles_poincare_surface(total_energy, time_integration, noOfInitialPoints);

  g1
      << "unset border\n"
      << "unset xtics\n"
      << "unset ytics\n"
      << "unset key\n";

  g1 << "set title \" Energy = 1/"
     << std::setprecision(1) << std::fixed
     << energy_denominator
     << "\"\n";
  g1<<std::setprecision(std::numeric_limits<double>::digits10 +1);
  g1 << "plot '-' with points pt 5 ps 0.001 lc rgb \"black\" \n";

  g1.send1d(std::make_tuple(y, py));
}

int main ()
{

  std::vector<float> energy_denominators{12, 11.2, 10.5, 10, 9.5, 8.8, 8, 7};

  const auto noOfInitialPoints = 3;
  const double time_integration = 400.0;
  Gnuplot g1;

  g1 << "set term pdfcairo enhanced font \"Arial,30\" size 11.7in,16.5in\n"
        "set out 'Henon_Helies2.pdf'\n";

  g1 << "set multiplot layout 4, 2 title \"Henon Heiles\" font \",45\n"
     << "set bmargin 5\n"
     << "set tmargin 2\n";

  for (const auto& energy_denominator : energy_denominators)
    plot_henon_heiles_poincare_surface(g1, energy_denominator, noOfInitialPoints, time_integration);

  g1 << "unset multiplot\n";

  return 0;

}