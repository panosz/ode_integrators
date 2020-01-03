from fractions import Fraction
from functools import total_ordering
import numpy as np


@total_ordering
class Resonance(Fraction):

    """Class for calculating the rotated dynamic system quantities
       near a given resonance ratio

       Resonance(ratio)

       Params:
       ratio: iterable or Fraction - like
            iterable must be equivalent (numerator, denominator)
       """

    def __new__(cls, numerator, denominator=None):
        parent_new = super(Resonance, cls).__new__
        try:
            self = parent_new(cls, numerator, denominator)
        except ValueError:
            msg = 'Invalid literal for {}: {}'
            raise ValueError(msg.format('Resonance', numerator))
        except TypeError:
            try:
                self = parent_new(cls, *numerator)
            except TypeError:
                msg = 'Invalid resonance ratio for {}: {}'
                raise TypeError(msg.format('Resonance', numerator))
        return self

    @property
    def ratio(self):
        return Fraction(self.numerator, self.denominator)

    def __repr__(self):
        cls = self.__class__.__name__
        return f'{cls}({self.ratio.numerator}, {self.ratio.denominator})'

    def jacobian(self):
        r = self.numerator
        s = self.denominator
        return np.array([[r, 0],
                         [-s, 1]])

    def transform_hessian(self, hessian):
        jac = self.jacobian()
        return jac.transpose() @ hessian @ jac

    def effective_mass(self, hessian):
        tr_hess = self.transform_hessian(hessian)
        return np.abs(tr_hess[0, 0])
