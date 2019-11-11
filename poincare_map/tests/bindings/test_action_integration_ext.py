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

        mass: scalar, positive: the mass of the oscillator.

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


class AnalyticPendulum(object):

    """Class with analytic expressions for the extended AnalyticPendulu"""

    def __init__(self, mass):
        """constructor for the analytic pendulum

        Parameters
        ----------
        mass: scalar, positive: the mass of the oscillator.

        """
        self._mass = mass

    def value(self, s):
        p = s[_P_COORDINATE]
        q = s[_Q_COORDINATE]
        F = s[_F_COORDINATE]

        return 0.5 * self._mass * p**2 - F * np.cos(q)


# Our algorithms assume a dynamics that is periodic in q.
# But Harmonic oscillator is not. We should therefore carefully
# Choose our initial points so that the span of the orbits
# is never too wide to reveal this non - periodicity
list_of_HO_points = [np.array(x) for x in [[10.1, 0.2, 1, 0],
                                           [3.1, 0.1, 2, 0],
                                           [2.1, 1.2, 3, 0],
                                           [1.1, 0.3, 4, 0],
                                           [2.1, 1.4, 5, 0],
                                           [2.1, 0.5, 6, 0],
                                           [2.1, 2.6, 7, 0],
                                           [2.1, 1.7, 8, 0],
                                           [2.1, 0.8, 9, 0],
                                           [1.1, -0.9, 10, 0],
                                           ]]
HO_mass = 1.0


@pytest.mark.parametrize("s", list_of_HO_points)
def test_action_integration(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'dt': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = HO_mass
    anal_ho = AnalyticHO(mass)

    ho_dynamic_system = ai.HarmonicOscDynamicSystem(mass)
    ho_action_integrals = ho_dynamic_system.action_integrals(s=s,
                                                             time=1000,
                                                             options=options)

    result = ho_action_integrals.hessian()
    desired = anal_ho.hessian(s)

    nt.assert_allclose(result, desired, atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("s", list_of_HO_points)
def test_orbit_closes(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'dt': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = HO_mass

    dynamic_system = ai.PendulumDynamicSystem(mass)

    a = dynamic_system.closed_orbit(s=s,
                                    time=1000,
                                    options=options)

    assert a.shape[1] == 4

    nt.assert_allclose(s[[0, 2]], a[-1, [0, 2]], atol=1e-10, rtol=1e-10)


@pytest.mark.parametrize("s", list_of_HO_points)
def test_following_orbit(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'dt': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0
    integration_time = 10

    dynamic_system = ai.PendulumDynamicSystem(mass)

    orbit, t = dynamic_system.orbit(s=s,
                                    time=integration_time,
                                    options=options)

    assert orbit.shape[1] == 4

    nt.assert_allclose(actual=t[-1], desired=integration_time)


@pytest.mark.parametrize("s", list_of_HO_points)
def test_orbit_energy_is_const(s):
    option_dict = {'abs_err': 1e-12, 'rel_err': 1e-12, 'dt': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    mass = 1.0

    dynamic_system = ai.PendulumDynamicSystem(mass)

    hamiltonian = dynamic_system.hamiltonian

    initial_energy = hamiltonian.value(s)

    a = dynamic_system.closed_orbit(s=s,
                                    time=1000,
                                    options=options)

    orbit_energies = np.apply_along_axis(lambda x: hamiltonian.value(x),
                                         axis=1,
                                         arr=a)

    nt.assert_allclose(actual=orbit_energies, desired=initial_energy)


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    def wrap_2pi(x):
        """
        maps number periodically on the range [0,2*pi)
        :param x: scalar or type convertible to numpy.array
        :return: scalar or numpy.array
        eg.

        >>> wrap_2pi([0, np.pi, 2 * np.pi])
        array([0.        , 3.14159265, 0.        ])

        >>> wrap_2pi(7)
        0.7168146928204138

        >>> wrap_2pi([3,4,5,6,7,8,9])
        array([3.        , 4.        , 5.        , 6.        , 0.71681469,
           1.71681469, 2.71681469])

        """
        return x - 2 * np.pi * np.floor_divide(x, 2 * np.pi)

    for p in list_of_HO_points:
        print(f"Energy of point{p} is {AnalyticPendulum(HO_mass).value(p)}")

    option_dict = {'abs_err': 1e-13, 'rel_err': 1e-13, 'dt': 1e-2}
    options = ai.IntegrationOptions(**option_dict)
    bad_point = list_of_HO_points[6]
    mass = 1.0
    dynamic_system = ai.PendulumDynamicSystem(mass)
    hamiltonian = dynamic_system.hamiltonian
    integration_time = 100

    print(f'bad_point = {bad_point}')
    print(dynamic_system.hamiltonian.value(bad_point))

    orbit, t = dynamic_system.orbit(s=bad_point,
                                    time=integration_time,
                                    options=options)

    orbit_energies = np.apply_along_axis(lambda x: hamiltonian.value(x),
                                         axis=1,
                                         arr=orbit)
    plt.plot(wrap_2pi(orbit[:, 1]), orbit[:, 0], 'k', markersize=.1)
    plt.show()
