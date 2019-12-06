import numbers
import numpy as np
import scipy as sp


def amplitudes(y):
    n_samples = y.shape[0]
    yf = sp.fft(y, axis=0)
    return yf / n_samples


class HarmonicSeries(object):
    """
    Custom container for Amplitudes of Harmonic Series.
    Attributes:
        data = an iterable (np.array) holding the amplidudes
        in standard DFT sequence

    Provides safe indexing for harmonic numbers.
    Raises exception if out of bounds negative or
    positive harmonic numbers are given as input.
    """
    def __init__(self, data):
        self.data = np.array(data)
        self.max_order = self.data.size // 2
        self.min_order = -((self.data.size - 1) // 2)

    def __len__(self):
        return self.data.size

    def __getitem__(self, key):
        return self.data[key]

    def __iter__(self):
        return iter(self.data)

    def _in_range(self, order):
        return order >= self.min_order and order <= self.max_order

    def nth_amplitude(self, order):
        if isinstance(order, numbers.Integral):
            if not self._in_range(order):
                msg = f'Harmonic amplitude of order {order} out of range'
                raise IndexError(msg)
            return self.data[order]
        else:
            order_type = order.__class__.__name__
            msg = f'Type {order_type} cannot be used as harmonic order.'
            raise TypeError(msg)


#  def sample_closed_orbit(closed_orbit, n_samples):
    #  t = np.linspace(0.0, closed_orbit.theta_period, n_samples, endpoint=False)
    #  return closed_orbit(t)


def orbital_harmonics_on_const_chi_bar(perturb, closed_orbit, n_samples):
    """
    samples a function along the projection of a closed_orbit
    at chi_bar = const.
    :param perturb: a callable with signature perturb(p,q,chi)
    :param closed_orbit: ds.Orbit type
    """

    delta_chi_bar = 1
    m = perturb.m

    p, q, f_chi, _ = sample_closed_orbit(closed_orbit, n_samples=n_samples)

    # prepare data for broadcasting
    p = p.reshape((n_samples, 1))
    q = q.reshape((n_samples, 1))
    f_chi = f_chi.reshape((n_samples, 1))

    chi_bar = np.array([0, delta_chi_bar]).reshape((1, 2))

    # sample perturb along the projection of the orbit on two different
    # chi_bar = const surfaces
    orbital_perturb = perturb(p=p, q=q, chi=f_chi + chi_bar)

    C_ns = amplitudes(orbital_perturb)
    A_ns = (C_ns[:, 1] * np.exp(1.0j * m) - C_ns[:, 0]) \
        / (np.exp(2.0j * m) - 1)

    return HarmonicSeries(A_ns)
