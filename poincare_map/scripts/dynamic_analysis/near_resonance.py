from fractions import Fraction
import numpy as np


class Resonance(object):

    """Class for calculating the rotated dynamic system quantities
       near a given resonance ratio

       Resonance(ratio)

       Params:
       ratio: iterable or Fraction - like
            iterable must be equivalent (numerator, denominator)
       """

    def __init__(self, numerator, denominator=None):
        try:
            self._ratio = Fraction(numerator, denominator)
        except ValueError:
            cls = type(self)
            msg = 'Invalid literal for {}: {}'
            raise ValueError(msg.format(cls, numerator))
        except TypeError:
            try:
                self._ratio = Fraction(*numerator)
            except TypeError:
                cls = type(self)
                msg = 'Invalid resonance ratio for {}: {}'
                raise TypeError(msg.format(cls, numerator))

    def __eq__(self, other):
        return self.ratio == other

    def __repr__(self):
        cls = self.__class__.__name__
        return f'{cls}({self.ratio.numerator}, {self.ratio.denominator})'

    def __str__(self):
        return str(self.ratio)

    @property
    def ratio(self):
        return self._ratio

    def jacobian(self):
        r = self._ratio.numerator
        s = self._ratio.denominator
        return np.array([[r, 0],
                         [-s, 1]])

    def transform_hessian(self, hessian):
        jac = self.jacobian()
        return jac.transpose() @ hessian @ jac

    def effective_mass(self, hessian):
        tr_hess = self.transform_hessian(hessian)
        return np.abs(tr_hess[0, 0])
