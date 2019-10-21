import numpy as np
import numpy.testing as nt
import pytest
import action_integration_ext as ai

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
                                        [4.1, 3.8 - np.pi, 9, 0],
                                        [4.1, 3.9 - np.pi, 10, 0],
                                        ]]


@pytest.mark.parametrize("s", list_of_points)
def test_action_integration(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'init_time_step': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0
    anal_ho = AnalyticHO(mass)
    numer_ho = ai.integrate_E_H_O(s,
                                  mass=mass,
                                  integration_time=1000,
                                  integration_options=options)
    result = numer_ho.hessian()
    desired = anal_ho.hessian(s)

    nt.assert_allclose(result, desired, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("s", list_of_points)
def test_action_integration2(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'init_time_step': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0
    anal_ho = AnalyticHO(mass)

    ho_dynamic_system = ai.HarmonicOscDynamicSystem(mass)
    ho_action_integrals = ho_dynamic_system.action_integrals(s=s,
                                                             integration_time=1000,
                                                             integration_options=options)

    result = ho_action_integrals.hessian()
    desired = anal_ho.hessian(s)

    nt.assert_allclose(result, desired, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("s", list_of_points)
def test_orbit_integration(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'init_time_step': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0

    a = ai.closed_pendulum_orbit(s,
                                 mass=mass,
                                 integration_time=1000,
                                 integration_options=options)

    assert a.shape[1] == 4

    nt.assert_allclose(s[[0, 2]], a[-1, [0, 2]], atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("s", list_of_points)
def test_orbit_integration2(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'init_time_step': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0

    ho_dynamic_system = ai.PendulumDynamicSystem(mass)

    a = ho_dynamic_system.closed_orbit(s=s,
                                       integration_time=1000,
                                       integration_options=options)

    assert a.shape[1] == 4

    nt.assert_allclose(s[[0, 2]], a[-1, [0, 2]], atol=1e-10, rtol=1e-10)
