//
// Created by Panagiotis Zestanakis on 01/03/19.
//

#ifndef ODE_INTEGRATORS_STATE_WRAPPER_HPP
#define ODE_INTEGRATORS_STATE_WRAPPER_HPP

#include <array>
#include <boost/range/algorithm/copy.hpp>
#include <boost/operators.hpp>
#include <boost/range/algorithm/max_element.hpp>

#include <boost/numeric/odeint/algebra/vector_space_algebra.hpp>

namespace DS
{
    class StateWrapper : boost::additive<StateWrapper, boost::additive<StateWrapper, double,
        boost::multiplicative<StateWrapper, double>
    > > {
      using T=std::array<double, 4>;
     private:
      T s_{0, 0, 0, 0};

     public:
      using value_type = typename T::value_type;
      using iterator =  typename T::iterator;
      using const_iterator = typename T::const_iterator;

      StateWrapper (value_type p, value_type q, value_type F, value_type chi)
          : s_{p, q, F, chi}
      { };
      StateWrapper ()
      { };

      inline value_type& operator[] (unsigned i)
      { return s_[i]; };

      inline value_type operator[] (unsigned i) const
      { return s_[i]; };

      inline auto begin ()
      {
        return std::begin(s_);
      }

      inline auto begin () const
      {
        return std::cbegin(s_);
      }

      inline auto end ()
      {
        return std::end(s_);
      }

      inline auto end () const
      {
        return std::cend(s_);
      }

      StateWrapper& operator+= (double d)
      {
        for (auto& x: s_)
          x += d;
        return *this;
      }

      StateWrapper& operator*= (double d)
      {
        for (auto& x: s_)
          x *= d;
        return *this;
      }

      StateWrapper& operator+= (const StateWrapper& other);

      StateWrapper operator/ (const StateWrapper& other) const;

      StateWrapper abs () const;

      double norm_inf () const;

    };

    StateWrapper abs (const StateWrapper& s);

    std::ostream& operator<< (std::ostream& os, const StateWrapper& state);

}

namespace boost { namespace numeric { namespace odeint {
            template<>
            struct vector_space_norm_inf< DS::StateWrapper >
            {
                typedef double result_type;
                double operator()( const DS::StateWrapper &s ) const
                {
                  return s.norm_inf();
                }
            };
        } } }

#endif //ODE_INTEGRATORS_STATE_WRAPPER_HPP
