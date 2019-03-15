//
// Created by Panagiotis Zestanakis on 06/03/19.
//

#include "input_output/hdf5_io.hpp"

arma::mat matrix_from_collection_of_armadillo_states (const std::vector<DS::armadillo_state>& states)
{
  arma::mat out;
  const auto number_of_states = states.size();

  if (number_of_states==0)
    return out;

  const auto number_of_dimensions = states[0].size();

  out.set_size(number_of_states, number_of_dimensions);

  for (long long unsigned i = 0; i < number_of_states; ++i)
    {
      out.row(i) = states[i];
    }
  return out;
}

std::ostream& operator<< (std::ostream& os, const std::vector<DS::armadillo_state>& states)
{
  {

    for (const auto& s : states)
      s.raw_print(os);
    return os;
  }
}

ArmaOrbitCrossOutput::ArmaOrbitCrossOutput (const OrbitCrossOutput<DS::armadillo_state>& orbitCrossOutput)
    :
    initial_point{orbitCrossOutput.initial_point},
    cross_points{matrix_from_collection_of_armadillo_states(orbitCrossOutput.cross_points)}
{}






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