//
// Created by Panagiotis Zestanakis on 05/03/19.
//

#include <boost/numeric/odeint.hpp>
#include <armadillo>
#include <cmath>
#include <boost/range/algorithm_ext/push_back.hpp>

using namespace arma;
using namespace boost::numeric::odeint;

using stateType = rowvec2;

using error_stepper_type = runge_kutta_cash_karp54<stateType,double, stateType, double, vector_space_algebra>;
using controlled_stepper_type = controlled_runge_kutta<error_stepper_type>;

namespace boost { namespace numeric { namespace odeint {
            template<>
            struct vector_space_norm_inf< stateType >
            {
                typedef double result_type;
                double operator()( const stateType &p ) const
                {
                  return norm(p,"inf");
                }
            };
        } } }



void h_o (const stateType& x, stateType& dxdt, double /*t*/)
{
  dxdt[0] = -sin(x[1]);
  dxdt[1] = x[0];
}

arma::mat fill_matrix_by_row (std::vector<stateType> rows)
{
  arma::mat out;

  if (!rows.empty())
  out.set_size(rows.size(),stateType::n_elem);

  for (int i = 0; i < rows.size(); ++i)
    {
      out.row(i) = rows[i];
    }
  return out;
}

int main ()
{

  stateType init = {0, 1};
  stateType derivs;
  std::cout << init;
  h_o(init,derivs,0);
  std::cout<<"derivative:\n"<<abs(derivs);

  const auto controlled_stepper = make_controlled(1e-8, 1e-8, error_stepper_type());

  auto orbit_iterators = make_adaptive_range(controlled_stepper,
                                             h_o,
                                             init, 0, 10, 1e-5);

  auto orbit_points = boost::make_iterator_range(orbit_iterators.first, orbit_iterators.second);


  std::cout<<"integrate orbit\n";

  std::vector<stateType> points;
  boost::push_back(points,orbit_points);

  mat Out = fill_matrix_by_row(points);

  std::cout<<Out;
}