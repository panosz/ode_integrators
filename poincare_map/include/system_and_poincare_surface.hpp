//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_SYSTEM_AND_POINCARE_SURFACE_HPP
#define ODE_INTEGRATORS_SYSTEM_AND_POINCARE_SURFACE_HPP

#include "myUtilities/wrap.hpp"
#include "phase_space_description.hpp"

struct Surface {
    CoordinateTag coordinate_tag;
    double position;
    int direction;
    Surface (CoordinateTag ct, double pos, int direct)
        : coordinate_tag(ct), position(pos), direction(direct)
    { }

    template <typename State>
    double eval(const State& s) const noexcept
    {
      const auto index = static_cast<unsigned>(coordinate_tag);
      const auto absolute_distance =  s[index] - position;

      if (coordinate_tag == CoordinateTag::p || coordinate_tag==CoordinateTag::F)
        return absolute_distance;

      else return PanosUtilities::wrap_minus_pi_pi(absolute_distance);

    }
    template <typename State>
    double operator()(const State& s) const noexcept
    {
      return eval(s);
    }

};

template<typename System>
class SystemAndPoincareSurface {
 private:
  System sys_;
  Surface surf_;

 public:
  using StateType = typename System::StateType;
  SystemAndPoincareSurface (const System& sys, Surface surf)
      : sys_{sys}, surf_{surf}
  { }


  void operator() (const StateType & s, StateType& dsdt, const double t) const
  {
    sys_(s, dsdt, t);
  }

  double surface_eval (const StateType& s) const
  {
    return surf_.eval(s);
  }

  Surface poincare_surface () const
  {
    return surf_;
  }

  void eval_perp_to_surface (const StateType& s, StateType& dsdt, const double t) const
  {
    sys_(s, dsdt, t);

    const auto index = static_cast<unsigned >(surf_.coordinate_tag);

    const auto normalization = dsdt[index];

    for (auto& item : dsdt)
      item = item / normalization;

  }
};

template <typename System>
SystemAndPoincareSurface<System> make_system_and_poincare_surface(System sys, Surface surf)
{
  return SystemAndPoincareSurface<System>{sys,surf};
}

#endif //ODE_INTEGRATORS_SYSTEM_AND_POINCARE_SURFACE_HPP
