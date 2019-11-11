#ifndef ORBIT_CLOSING_CHECK_HPP_WWU3ENM5
#define ORBIT_CLOSING_CHECK_HPP_WWU3ENM5
#include "integration_options.hpp"
#include "phase_space_description.hpp"

namespace OrbitClosing
{


  /// return double 'half_width'. half_width is the actual tolerance
  /// around 'y', for which a float 'x' should compare almost
  /// equal to 'y'
  /// half_width = 'atol' + 'rtol' * abs('y')
  double tolerance_half_width(double y,
      double atol,
      double rtol) noexcept;

  /// return true if 'x' and 'y' are within the specified tolerance
  /// abs('x' - 'y') <= 'atol' + 'rtol' * abs('y')
  bool almost_equal(double x,
      double y,
      double atol,
      double rtol) noexcept;


  /// return true if 'x' and 'y' are within the specified tolerance
  /// abs('x' - 'y') <= 'half_width'
  bool almost_equal(double x,
      double y,
      double half_width ) noexcept;


  class P_Checker
  {
    /// callable. When called as occ(point), returns true,
    /// when the p coordinate of 'point' is close to the
    /// p coordinate of 'p_init', the point that was used
    /// to initialize occ.

    private:
      double p_{};
      double half_width_{};

    public:
      P_Checker(double init_p, double half_width);

      template <typename PointType>
        bool operator()(const PointType& point) const noexcept
        {
          const auto p = point[static_cast<unsigned>(CoordinateTag::p)];
          return almost_equal(p, p_, half_width_);
        }
  };

  template<typename PointType>
    P_Checker
    make_P_Checker(const PointType& init_point, const IntegrationOptions& options)
    {
      const auto init_p = init_point[static_cast<unsigned>(CoordinateTag::p)];
      const auto & coef = options.coefficient;
      const auto & atol = options.abs_err;
      const auto & rtol = options.rel_err;

      const auto half_width = coef * tolerance_half_width(init_p, atol, rtol);

      return P_Checker(init_p, half_width);
    }

}




#endif /* end of include guard: ORBIT_CLOSING_CHECK_HPP_WWU3ENM5 */
