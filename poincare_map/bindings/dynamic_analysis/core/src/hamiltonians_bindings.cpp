#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "hamiltonians_bindings.hpp"
#include "UnperturbedExtendedOscillatorHamiltonian.hpp"
#include "UnperturbedExtendedPendulumHamiltonian.hpp"


namespace p = boost::python;
namespace np = boost::python::numpy;


namespace HamiltoniansBindings{



  template<typename Hamiltonian>
    void export_hamiltonian(const char* hamiltonian_name,
        const char* docstring="")
    {
      using HD = Decorator<Hamiltonian>;

      p::class_<HD>(hamiltonian_name,
          docstring,
          p::init<HD>())
        .def("value",
            &HD::value,
            (p::arg("s")));
    }


  void export_hamiltonians()
  {
    using HarmOscHamiltonian = DS::UnperturbedExtendedOscillatorHamiltonian;
    using PendulumHamiltonian = DS::UnperturbedExtendedPendulumHamiltonian;
    export_hamiltonian<HarmOscHamiltonian>("H_O_Hamiltonian");
    export_hamiltonian<PendulumHamiltonian>("PendulumHamiltonian");
  }

}

