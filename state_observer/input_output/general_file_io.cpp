//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#include <iostream>
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



auto make_file_paths(const char *filename, const char* extension)
{
  namespace fs = std::filesystem;

  auto path = fs::path(filename);
  path += extension;

  return path;
}

auto make_text_file_paths(const char *filename)
{
  return make_file_paths(filename,".txt");
}

auto make_hdf5_file_paths(const char *filename)
{
  return make_file_paths(filename,".hdf5");
}

std::ofstream prepare_text_files_for_output (const char *filename)
{
  auto path = make_text_file_paths(filename);

  std::ofstream outputFile(path.native());
  if (!outputFile)
    throw std::runtime_error("cannot open output file for exact points");


  outputFile << std::setprecision(16);

  return outputFile;
}

void remove_existing_hdf5_file(const std::filesystem::path& filepath)
{
  namespace fs = std::filesystem;

  if(fs::exists(filepath))
    {
      std::cerr<<"Warning replacing file "<< filepath<<std::endl;
      fs::remove(filepath);
    }
}

std::filesystem::path prepare_hdf5_files_for_output(const char * filename)
{

  namespace fs = std::filesystem;
  auto path = make_hdf5_file_paths(filename);

  remove_existing_hdf5_file(path);

  return path;
}