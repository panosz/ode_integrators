#!/usr/bin/env python3
# Using the doctest module here to ensure that the results are as expected.
r'''>>> from hello_ext import *
    >>> args_to_array(1.,2.,3.)
    array([1., 2., 3.])
    >>> a = array_builder()
    >>> a.to_array(1,2,3)
    array([1., 2., 3.])
    >>> args_to_Hessian(1,2,3)
    array([[1., 2.],
           [2., 3.]])
    >>> a = make_myPoint(3,4)
    >>> type(a)
    <class 'hello_ext.MyPythonPoint'>
    >>> a.x
    3.0
    >>> a.x=3.14
    >>> a.x
    3.14
    >>> a.y
    4.0
'''


if __name__ == '__main__':
    import doctest
    doctest.testmod()
