//
// Created by Panagiotis Zestanakis on 24/05/18.
//

#ifndef ODE_INTEGRATORS_UTILITIES_HPP
#define ODE_INTEGRATORS_UTILITIES_HPP

namespace PanosOde::Utils
{

    template <typename State>
    class PoincareSurfaceCrossEventObserver {
     protected:
      mutable bool first_call = true;
      mutable double previous_distance_from_surface_ = 0;

      virtual double distance_from_surface (const State& ) const =0;
      virtual bool event_predicate (double ) const = 0; // distance decreases, only when we cross the surface

      void update_distance (double distance) const
      {
        previous_distance_from_surface_ = distance;
      }

     public:
      bool operator() (const State& x, double /*t*/) const
      {
        if (!first_call)
          {
            const double distance = distance_from_surface(x);
            const bool crossed_surface = event_predicate(distance);
            update_distance(distance);
            return crossed_surface;
          }
        else
          {
            first_call = false;
            update_distance(distance_from_surface(x));
            return false;
          }

      }

      double distance () const
      {
        return previous_distance_from_surface_;
      }

    };




    template<typename StateFunctor, typename EventObserver, typename Stepper, typename System>
    class EventFunctor {
     private:
      const StateFunctor sf_;
      const EventObserver eo_;
      Stepper stepper_;
      System sys_;

     public:
      EventFunctor (StateFunctor sf, EventObserver eo, Stepper stepper, System sys)
          :
          sf_{sf}, eo_{eo}, stepper_{stepper}, sys_{sys}
      { };

      template<typename State>
      void operator() (const State& x, double t)
      {
        if (eo_(x, t))
          {
            double distance_from_surf = eo_.distance();
            State on_surface_state;
            State derivs;
            sys_(x, derivs, t);
            stepper_.do_step(sys_, x, derivs, t, on_surface_state, -distance_from_surf);
            sf_(on_surface_state, t);
//        sf_(x, t);
          }
      }

    };

    template<typename StateFunctor, typename EventObserver, typename Stepper, typename System>
    EventFunctor<StateFunctor, EventObserver, Stepper, System>
    makeEventFunctor (StateFunctor sf, EventObserver eo, Stepper st, System sy)
    {
      return EventFunctor<StateFunctor, EventObserver, Stepper, System>(sf, eo, st, sy);
    }



}

#endif //ODE_INTEGRATORS_UTILITIES_HPP
