//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#include "input_output/text_file_io.hpp"

void write_to_text_files (const char *filename, const OrbitCrossOutput& crossOutput)
{

  auto outputfile = prepare_text_files_for_output(filename);

  outputfile << crossOutput;

  outputfile.close();

}
void write_to_text_files (const char *filename, const std::vector<OrbitCrossOutput>& crossOutputs)
{
  auto exact_file = prepare_text_files_for_output(filename);

  for (const auto& out : crossOutputs)
    {
      exact_file << out << '\n';
    }

  exact_file.close();

}
