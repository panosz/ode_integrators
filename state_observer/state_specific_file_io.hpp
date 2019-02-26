//
// Created by Panagiotis Zestanakis on 26/02/19.
//
#ifndef ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
#define ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP

#include "orbit_classes.hpp"
#include "general_file_io.hpp"

template <typename State>
void write_to_files (const char *filename, const ApproximateAndExactCrossOutput<State>& crossOutput)
{

  auto[exact_file, approx_file] = prepare_files_for_output(filename);

  exact_file << crossOutput.exact;
  approx_file << crossOutput.approximate;

  exact_file.close();
  approx_file.close();

}

template <typename State>
void write_to_files (const char *filename, const std::vector<ApproximateAndExactCrossOutput<State>>& crossOutputs)
{
  auto[exact_file, approx_file] = prepare_files_for_output(filename);

  for (const auto& out : crossOutputs)
    {
      exact_file << out.exact << '\n';
      approx_file << out.approximate << '\n';
    }

  exact_file.close();
  approx_file.close();

}

#endif //ODE_INTEGRATORS_STATE_SPECIFIC_FILE_IO_HPP
