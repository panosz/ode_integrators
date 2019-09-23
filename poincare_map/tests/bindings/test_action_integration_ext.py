import numpy as np
import numpy.testing as nt
import pytest
import action_integration_ext as ai


def test_args_to_array():
    result = ai.args_to_array(1., 2., 3.)
    desired = np.array([1., 2., 3.])
    nt.assert_array_equal(result, desired)


def test_args_to_Hessian():
    result = ai.args_to_Hessian(1, 2, 3)
    desired = np.array([[1., 2.],
                        [2., 3.]])
    nt.assert_array_equal(result, desired)


_P_COORDINATE = 0
_Q_COORDINATE = 1
_F_COORDINATE = 2
_PHI_COORDINATE = 3


class AnalyticHO(object):

    """Class with analytic expressions for the extended Harmonic Oscillator
    in Action Space"""

    def __init__(self, mass):
        """construct the harmonic oscillator

        :mass: scalar, positive: the mass of the oscillator.

        """
        self._mass = mass

    def value(self, s):
        p = s[_P_COORDINATE]
        q = s[_Q_COORDINATE]
        F = s[_F_COORDINATE]

        return self._mass * p**2 + F * q**2

    def action(self, s):
        F = s[_F_COORDINATE]
        return self.value(s) / (2 * np.sqrt(self._mass * F))

    def dKdJ(self, s):
        F = s[_F_COORDINATE]
        return 2 * np.sqrt(self._mass * F)
        pass

    def dKdF(self, s):
        F = s[_F_COORDINATE]
        return self.action(s) * np.sqrt(self._mass / F)

    def d2KdJ2(self, s):
        return 0

    def d2KdJdF(self, s):
        F = s[_F_COORDINATE]
        return np.sqrt(self._mass / F)

    def d2KdF2(self, s):
        F = s[_F_COORDINATE]
        J = self.action(s)
        return -J * np.sqrt(self._mass / F) / (2 * F)

    def hessian(self, s):

        return np.array([[self.d2KdJ2(s), self.d2KdJdF(s)],
                         [self.d2KdJdF(s), self.d2KdF2(s)]])


list_of_points = [np.array(x) for x in [[4.1, 3.0, 1, 0],
                                        [4.1, 3.1, 2, 0],
                                        [4.1, 3.2, 3, 0],
                                        [4.1, 3.3, 4, 0],
                                        [4.1, 3.4, 5, 0],
                                        [4.1, 3.5, 6, 0],
                                        [4.1, 3.6, 7, 0],
                                        [4.1, 3.7, 8, 0],
                                        [4.1, 3.8, 9, 0],
                                        [4.1, 3.9, 10, 0],
                                        ]]


@pytest.mark.parametrize("s", list_of_points)
def test_action_integration(s):
    mass = 1.0
    anal_ho = AnalyticHO(mass)
    numer_ho = ai.integrate_E_H_O(s, mass=mass)
    result = numer_ho.hessian()
    desired = anal_ho.hessian(s)
    nt.assert_allclose(result, desired, atol=1e-10, rtol=1e-10)
