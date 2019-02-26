//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_SYSTEM_AND_POINCARE_SURFACE_HPP
#define ODE_INTEGRATORS_SYSTEM_AND_POINCARE_SURFACE_HPP


enum class VariableTag : unsigned {
  p = 0,
  q = 1,
  chi = 2
};

struct Surface {
    VariableTag variable_tag;
    double position;
    int direction;
    Surface (VariableTag vt, double pos, int direct)
        : variable_tag(vt), position(pos), direction(direct)
    { }
};

template<typename System>
class SystemAndPoincareSurface {
 private:
  System sys_;
  const Surface surf_;

 public:
  using StateType = typename System::StateType;
  SystemAndPoincareSurface (const System& sys, Surface surf)
      : sys_{sys}, surf_{surf}
  { }

  void eval (const StateType & s, StateType& dsdt, const double t) const
  {
    sys_(s, dsdt, t);
  }

  double surface_eval (const StateType& s) const
  {
    const auto index = static_cast<unsigned>(surf_.variable_tag);

    return s[index] - surf_.position;
  }

  Surface poincare_surface () const
  {
    return surf_;
  }

  void eval_perp_to_surface (const StateType& s, StateType& dsdt, const double t) const
  {
    sys_(s, dsdt, t);

    const auto index = static_cast<unsigned >(surf_.variable_tag);

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
