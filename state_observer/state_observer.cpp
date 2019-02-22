//
// Created by Panagiotis Zestanakis on 18/09/18.
//

#include <iostream>
#include <array>
#include <vector>
#include <boost/numeric/odeint.hpp>
#include "../utilities.hpp"

using namespace boost::numeric::odeint;

/* The type of container used to hold the state vector */
typedef std::array<double, 2> State;

class System {
 public:
  System () = default;

  void operator() (const State& x, State& dxdt, const double /*t*/)
  {
    dxdt[0] = 1;
    dxdt[1] = x[0];

  }
};

class SystemCross
{
 public:
  SystemCross()=default;

  void operator() (const State& x, State& dxdt, const double /*t*/)
  {
    dxdt[0] = 1/x[0];
    dxdt[1] = 1;

  }

};

class Value_Cross_Event_Observer : public PanosOde::Utils::PoincareSurfaceCrossEventObserver<State> {
  const std::vector<double>& values;
  mutable std::vector<double>::const_iterator value_iter;

 public:
  Value_Cross_Event_Observer (const std::vector<double>& val)
      :values{val},value_iter{values.cbegin()}
  { };
 protected:
  double distance_from_surface (const State& state) const override
  {
    return state[1] - *value_iter;
  }
 private:
  bool event_predicate (double distance) const override
  {
    if(value_iter != values.cend())
      {
        const auto pred = (distance > 0) && (previous_distance_from_surface_ < 0);
        if (pred)
          {
            std::cout<<"crossed_surf at "<<*value_iter<<"\n";
            ++value_iter;
          }
        return pred;
      }
    return false;
  }

};

using ErrorStepperType = runge_kutta_cash_karp54<State>;
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType>;

int main ()
{
  State init_state{0, 0};

  std::vector<double> psi;

  const auto push_back_observer = [&psi] (State& state, double /*t*/)
  { psi.push_back(state[1]); };

  const std::vector<double> values = {1, 2,3, 4,5,6,7,8,9, 10};

  const double abs_err = 1.0e-16;
  const double rel_err = 1.0e-16;
  const double max_dt = 0.001;
  const auto controlled_stepper = make_controlled(abs_err, rel_err,max_dt, ErrorStepperType());

  auto cross_event_observer = PanosOde::Utils::makeEventFunctor(push_back_observer, Value_Cross_Event_Observer(values),
                                                                 ErrorStepperType(), SystemCross());

  integrate_adaptive(controlled_stepper, System(), init_state, 0.0, 10.0, 1e-5, cross_event_observer);


  for (auto out:psi)
    std::cout << out << '\n';

}
