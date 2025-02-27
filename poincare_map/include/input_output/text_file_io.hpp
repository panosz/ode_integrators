//
// Created by Panagiotis Zestanakis on 26/02/19.
//
#ifndef ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
#define ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP

#include "samplingCollections.hpp"
#include "prepare.hpp"
#include "myUtilities/data_reading.hpp"

template<typename State>
std::vector<State> get_state_from_file (const char *input_filename, unsigned number_of_dimensions)
{
  std::vector<State> states{};

  auto input_file = open_input_file(input_filename);


  for (std::string line; getline(input_file, line);)
    {
      auto cleared_line = PanosUtilities::trimm_comments(line, "#");
      if (!cleared_line.empty())
        {
          const auto state_coordinates = PanosUtilities::doubles_from_string(cleared_line);

          if (state_coordinates.size() != number_of_dimensions)
            throw std::runtime_error("get_state_from_file: Unexpected number of coordinates, expecting " +
                                     std::to_string(number_of_dimensions)+ " got "+
                                     std::to_string(state_coordinates.size()));

          State  input_state;
          boost::copy(state_coordinates,input_state.begin());
          states.push_back(input_state);
        }

    }

  return states;
}



void write_to_text_files (const char *filename, const OrbitCrossOutput& crossOutput);

void write_to_text_files (const char *filename, const std::vector<OrbitCrossOutput>& crossOutputs);

#endif //ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
