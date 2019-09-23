#!/usr/bin/env python3
# Using the doctest module here to ensure that the results are as expected.
r'''>>> from action_integration_ext import *
    >>> args_to_array(1.,2.,3.)
    array([1., 2., 3.])
    >>> args_to_Hessian(1,2,3)
    array([[1., 2.],
           [2., 3.]])
    >>> b=make_action_integration_result()
    >>> b.hessian()
    array([[-0.,  0.],
           [ 0.,  0.]])
'''


if __name__ == '__main__':
    import doctest
    doctest.testmod()
