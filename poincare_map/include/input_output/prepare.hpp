//
// Created by Panagiotis Zestanakis on 26/02/19.
//

#ifndef ODE_INTEGRATORS_FILE_IO_HPP
#define ODE_INTEGRATORS_FILE_IO_HPP

#include <filesystem>
#include <fstream>

std::ifstream open_input_file(const char *input_filename);



std::ofstream prepare_text_files_for_output (const char *filename);


std::filesystem::path prepare_hdf5_files_for_output(const char * filename);



#endif //ODE_INTEGRATORS_FILE_IO_HPP
