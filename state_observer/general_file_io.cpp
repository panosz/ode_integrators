//
// Created by Panagiotis Zestanakis on 26/02/19.
//


#include <iomanip>

#include "general_file_io.hpp"
std::ifstream open_input_file (const char *input_filename)
{


  namespace fs = std::filesystem;
  auto path = fs::path(input_filename);

  if(!fs::exists(path))
    {
      throw std::runtime_error(std::string("input file does not exist: ") + input_filename);
    }


  std::ifstream input_file(path.native());
  if (!input_file.is_open())
    throw std::runtime_error(std::string("cannot open input file") + input_filename);

  return input_file;

}

auto make_file_paths (const char *filename)
{
  namespace fs = std::filesystem;

  auto exact_path = fs::path(filename);
  auto approx_path = exact_path;

  exact_path += "_exact.txt";
  approx_path += "_approx.txt";

  return std::make_pair(exact_path, approx_path);
}


std::pair<std::ofstream, std::ofstream> prepare_files_for_output (const char *filename)
{
  auto[exact_path, approx_path] = make_file_paths(filename);

  std::ofstream exact_file(exact_path.native());
  if (!exact_file)
    throw std::runtime_error("cannot open output file for exact points");

  std::ofstream approx_file(approx_path.native());
  if (!exact_file)
    throw std::runtime_error("cannot open output file for approximate points");

  exact_file << std::setprecision(16);
  approx_file << std::setprecision(16);

  return std::make_pair(std::move(exact_file),
                        std::move(approx_file));//std::move required, because ofsream copy constructor is deleted
}

