//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#include "hdf5_io.hpp"

void write_as_group(const std::string& filename, const std::string& group, const ArmaOrbitCrossOutput& armaOrbitCrossOutput )
{
  using namespace arma;

  armaOrbitCrossOutput.initial_point.save(hdf5_name(filename,group+"/initial_point",hdf5_opts::append + hdf5_opts::trans));
  armaOrbitCrossOutput.cross_points.save(hdf5_name(filename, group+"/cross_points",hdf5_opts::append + hdf5_opts::trans));
}

void write_vector_of_datasets( std::filesystem::path output_path,
                               const std::string& group_prefix,
                               const std::vector<OrbitCrossOutput<DS::armadillo_state>>& crossOutputs)
{


  for (long unsigned int i = 0; i < crossOutputs.size(); ++i)
    {

      std::string group= group_prefix + std::to_string(i);
      std::cout<<"writing "<<group<<'\n';

      write_as_group(output_path,group, ArmaOrbitCrossOutput{crossOutputs[i]});

    }

}




void write_to_hdf5_files (const char *filename,
                          const std::vector<OrbitCrossOutput<DS::armadillo_state>>& crossOutputs)
{
  auto paths = prepare_hdf5_files_for_output(filename);

  write_vector_of_datasets(paths,"orbit",crossOutputs);

}