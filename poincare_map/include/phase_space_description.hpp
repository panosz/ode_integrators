#ifndef PHASE_SPACE_DESCRIPTION_HPP_YPKWLGNY
#define PHASE_SPACE_DESCRIPTION_HPP_YPKWLGNY


const unsigned PHASE_SPACE_VARIABLES = 4;
const unsigned PATH_INTEGRALS = 8;
const unsigned EXTENDED_SPACE_VARIABLES = PHASE_SPACE_VARIABLES + PATH_INTEGRALS;

enum class CoordinateTag : unsigned {
  p = 0,
  q = 1,
  F = 2,
  phi = 3,
  J = 4,
  t = 5,
  beta = 6,
  gamma = 7,
  beta1 = 8,
  gamma1 = 9,
  beta2 = 10,
  gamma2 =11
};

#endif /* end of include guard: PHASE_SPACE_DESCRIPTION_HPP_YPKWLGNY */
