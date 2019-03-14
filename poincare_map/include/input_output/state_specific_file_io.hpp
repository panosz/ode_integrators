//
// Created by Panagiotis Zestanakis on 26/02/19.
//
#ifndef ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
#define ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP

#include "samplingCollections.hpp"
#include "general_file_io.hpp"
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



template <typename State>
void write_to_text_files (const char *filename, const OrbitCrossOutput<State>& crossOutput)
{

  auto outputfile = prepare_text_files_for_output(filename);

  outputfile << crossOutput;

  outputfile.close();

}

template <typename State>
void write_to_text_files (const char *filename, const std::vector<OrbitCrossOutput<State>>& crossOutputs)
{
  auto exact_file = prepare_text_files_for_output(filename);

  for (const auto& out : crossOutputs)
    {
      exact_file << out << '\n';
    }

  exact_file.close();

}

#endif //ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
