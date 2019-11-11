//
//Created by Panagiotis Zestanakis on Nov 11, 2019
//
#ifndef UNPERTURBEDEXTENDEDOSCILLATORHAMILTONIAN_HPP_ZFQ3IXZ4
#define UNPERTURBEDEXTENDEDOSCILLATORHAMILTONIAN_HPP_ZFQ3IXZ4
#include "fields_and_brackets.hpp"
#include "phase_space_state_types.hpp"
namespace DS
{
    class UnperturbedExtendedOscillatorHamiltonian {

     private:
      double M_;

     public:

      explicit UnperturbedExtendedOscillatorHamiltonian (double M)
          : M_{M}
      { };

      template<typename ST>
      double value(const ST& s) const noexcept;

      template<typename ST>
      double operator() (const ST& s) const noexcept;

      template<typename ST>
      FirstDerivatives first_derivatives (const ST& s) const noexcept;

      template<typename ST>
      SecondDerivatives second_derivatives (const ST& s) const noexcept;

      template<typename ST>
      ThirdDerivatives third_derivatives (const ST& /*s*/) const noexcept;

      double action(const PhaseSpaceState& s) const;

      double dKdJ(const PhaseSpaceState& s) const noexcept;

      double dKdF(const PhaseSpaceState& s) const;

      double d2KdJ2(const PhaseSpaceState& /*s*/) const noexcept;

      double d2KdJdF(const PhaseSpaceState& s) const;

      double d2KdF2(const PhaseSpaceState& s) const;


    };
}

#endif /* end of include guard: UNPERTURBEDEXTENDEDOSCILLATORHAMILTONIAN_HPP_ZFQ3IXZ4 */
