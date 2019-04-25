//
// Created by Panagiotis Zestanakis on 23/04/19.
//

#include "hamiltonian_dynamic_system.hpp"

DS::FieldAndFirstDerivatives::FieldAndFirstDerivatives (const DS::Field& F, const DS::FieldFirstDerivatives& dF)
    : f{F}, df{dF}
{ }
