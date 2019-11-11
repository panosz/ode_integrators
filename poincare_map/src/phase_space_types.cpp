#include "phase_space_state_types.hpp"
namespace DS
{
  ExtendedSpaceState phase_to_extended_space_state (const PhaseSpaceState& pss)
  {
    ExtendedSpaceState ess{};
    ess.zeros();

    ess[static_cast<unsigned>(CoordinateTag::p)] =
      pss[static_cast<unsigned>(CoordinateTag::p)];
    ess[static_cast<unsigned>(CoordinateTag::q)] =
      pss[static_cast<unsigned>(CoordinateTag::q)];
    ess[static_cast<unsigned>(CoordinateTag::F)] =
      pss[static_cast<unsigned>(CoordinateTag::F)];
    ess[static_cast<unsigned>(CoordinateTag::phi)] =
      pss[static_cast<unsigned>(CoordinateTag::phi)];

    return ess;
  }
}
