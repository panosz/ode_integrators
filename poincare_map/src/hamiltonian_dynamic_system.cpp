//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#include "hamiltonian_dynamic_system.hpp"

namespace DS
{

   double UnperturbedExtendedOscillatorHamiltonian::action (const PhaseSpaceState& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return value(s)/(2 * std::sqrt(F));
    }

    double UnperturbedExtendedOscillatorHamiltonian::dKdJ (const PhaseSpaceState& s) const noexcept
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return 2 * std::sqrt(F);
    }

    double UnperturbedExtendedOscillatorHamiltonian::dKdF (const PhaseSpaceState& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return action(s) / std::sqrt(F);
    }

    double UnperturbedExtendedOscillatorHamiltonian::d2KdJ2 (const PhaseSpaceState&) const noexcept
    {
      return 0;
    }

    double UnperturbedExtendedOscillatorHamiltonian::d2KdJdF (const PhaseSpaceState& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      return 1 / std::sqrt(F);
    }


    double UnperturbedExtendedOscillatorHamiltonian::d2KdF2 (const PhaseSpaceState& s) const
    {
      const auto& F = s[static_cast<unsigned>(CoordinateTag::F)];

      const auto J = action(s);

      return -J / (2 * std::pow(F, 1.5));
    }


}