#ifndef PHASE_SPACE_STATE_TYPES_HPP_PZWRNE6F
#define PHASE_SPACE_STATE_TYPES_HPP_PZWRNE6F

#include "armadillo_state.hpp"
#include "phase_space_description.hpp"

namespace DS
{
    using PhaseSpaceState = armadillo_state<PHASE_SPACE_VARIABLES>;
    using ExtendedSpaceState = armadillo_state<EXTENDED_SPACE_VARIABLES>;

    ExtendedSpaceState phase_to_extended_space_state(const PhaseSpaceState &pss);
}


#endif /* end of include guard: PHASE_SPACE_STATE_TYPES_HPP_PZWRNE6F */
