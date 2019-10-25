#include <cmath>
#include "orbit_closing_check.hpp"

namespace OrbitClosing
{
    /// return double 'half_width'. half_width is the actual tolerance
    /// around 'y', for which a float 'x' should compare almost
    /// equal to 'y'
    /// half_width = 'atol' + 'rtol' * abs('y')
    double tolerance_half_width(double y,
                                double atol,
                                double rtol) noexcept
    {
      return atol + rtol * std::abs(y);
    }

    /// return true if 'x' and 'y' are within the specified tolerance
    /// abs('x' - 'y') <= 'atol' + 'rtol' * abs('y')
    bool almost_equal(double x,
                      double y,
                      double atol,
                      double rtol) noexcept
    {
      return std::abs(x - y) <= tolerance_half_width(y, atol, rtol) ;
    }


    /// return true if 'x' and 'y' are within the specified tolerance
    /// abs('x' - 'y') <= 'half_width'
    bool almost_equal(double x,
                      double y,
                      double half_width ) noexcept
    {
      return std::abs(x - y) <= half_width;
    }


    P_Checker::P_Checker(double init_p, double half_width)
            :p_{init_p},
            half_width_{half_width}
          {}
}
