//
// Created by Panagiotis Zestanakis on 01/03/19.
//

#include "state_wrapper.hpp"
DS::StateWrapper& DS::StateWrapper::operator+= (const DS::StateWrapper& other)
{
  for (unsigned long i = 0; i < std::size(s_); ++i)
    s_[i] += other.s_[i];
  return *this;
}
DS::StateWrapper DS::StateWrapper::operator/ (const DS::StateWrapper& other) const
{
  StateWrapper output{*this};
  for (unsigned long i = 0; i < std::size(s_); ++i)
    output.s_[i] /= other.s_[i];
  return output;
}
DS::StateWrapper DS::StateWrapper::abs () const
{
  StateWrapper output{*this};
  for (auto& x: output.s_)
    x = std::abs(x);
  return output;
}
double DS::StateWrapper::norm_inf () const
{
  const auto s_abs = abs();
  auto max_it = boost::range::max_element(s_abs.s_);
  return *max_it;
}
DS::StateWrapper DS::abs (const DS::StateWrapper& s)
{
  return s.abs();
}
std::ostream& DS::operator<< (std::ostream& os, const DS::StateWrapper& state)
{
  {
    boost::copy(state, std::ostream_iterator<typename StateWrapper::value_type>(os, " "));
    return os;
  }
}
