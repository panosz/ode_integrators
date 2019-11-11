//
//Created by Panagiotis Zestanakis on Nov 11, 2019
//
#ifndef UNPERTURBEDEXTENDEDPENDULUMHAMILTONIAN_HPP_LHYNB7IR
#define UNPERTURBEDEXTENDEDPENDULUMHAMILTONIAN_HPP_LHYNB7IR

#include "fields_and_brackets.hpp"
#include "phase_space_state_types.hpp"
namespace DS
{

    class UnperturbedExtendedPendulumHamiltonian {

     private:
      double M_;

     public:
      explicit UnperturbedExtendedPendulumHamiltonian (double M)
          : M_{M}
      { };

      template<typename ST>
      double value(const ST& s) const noexcept;

      template<typename ST>
      double operator() (const ST& s) const;

      template<typename ST>
      FirstDerivatives first_derivatives (const ST& s) const noexcept;

      template<typename ST>
      SecondDerivatives second_derivatives (const ST& s) const noexcept;

      template<typename ST>
      ThirdDerivatives third_derivatives (const ST& s) const noexcept;

    };
}

#endif /* end of include guard: UNPERTURBEDEXTENDEDPENDULUMHAMILTONIAN_HPP_LHYNB7IR */
