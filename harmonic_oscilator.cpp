#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/numeric/odeint.hpp>
#include "gnuplot-iostream.h"

using namespace boost::numeric::odeint;


/* The type of container used to hold the state vector */
typedef std::vector< double > State;

/* The rhs of x' = f(x) defined as a class */
class HarmonicOsc {

  double m_gam = 0;

 public:
  HarmonicOsc( double gam ) : m_gam(gam) { }

  void operator() ( const State &x , State &dxdt , const double /* t */ )
  {
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - m_gam*x[1];
  }
};

using BasicStepperType=runge_kutta4< State >;
using ErrorStepperType = runge_kutta_cash_karp54< State >;
using ControlledStepperType = controlled_runge_kutta<ErrorStepperType>;

class saving_observer
{
 private:
    std::vector<double> & xout_;
    std::vector<double> & tout_;
    const size_t every_=100;
    mutable size_t count = 0;
 public:
  saving_observer(std::vector<double>& xout, std::vector<double>& tout, size_t every):
      xout_{xout},tout_{tout},every_{every>0?every:1}
  {};
  void operator() (const State & x, double t) const
  {
    if(!(count++%every_))
      {
        xout_.push_back(x[0]);
        tout_.push_back(t);
      }

  };

};

int main ()
{


  Gnuplot g1;
  HarmonicOsc harmonic_oscillator(0.15);

  g1<<"plot '-' with lines, '-' with lines, '-' with lines \n";


  {//integrate with const step
    State x;
    x.emplace_back(0);
    x.emplace_back(1);
    BasicStepperType stepper;
    std::vector<double> t_out;
    std::vector<double> x_out;
    integrate_const( stepper , harmonic_oscillator , x , 0.0 , 100.0 , 0.1, saving_observer(x_out,t_out,10) );
    std::cout<< "constant integration output: "<<x[0]<< ", "<<x[1]<<std::endl;
    g1.send1d(std::make_tuple(t_out,x_out));
  }

  {//integrate with adaptive step using the explicitly declared ControlledStepperType
    State x;
    x.emplace_back(0);
    x.emplace_back(1);
    ControlledStepperType controlled_stepper;
    std::vector<double> t_out;
    std::vector<double> x_out;
    integrate_adaptive( controlled_stepper , harmonic_oscillator , x , 0.0 , 100.0 , 0.01,saving_observer(x_out,t_out,10) );
    std::cout<< "adaptive integration output(1): "<<x[0]<< ", "<<x[1]<<std::endl;
    g1.send1d(std::make_tuple(t_out,x_out));

  }

  {//integrate with adaptive step using make_controlled
    State x;
    x.emplace_back(0);
    x.emplace_back(1);
    double abs_err = 1.0e-16;
    double rel_err = 1.0e-14;
    auto controlled_stepper = make_controlled(abs_err,rel_err,ErrorStepperType());
    std::vector<double> t_out;
    std::vector<double> x_out;

    integrate_adaptive( controlled_stepper , harmonic_oscillator , x , 0.0 , 100.0 , 0.01,saving_observer(x_out,t_out,10) );
    std::cout<< "adaptive integration output(2): "<<x[0]<< ", "<<x[1]<<std::endl;
    g1.send1d(std::make_tuple(t_out,x_out));

  }

  int a;
  std::cin>> a;
  return 0;
}