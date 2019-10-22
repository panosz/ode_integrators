import numpy as np
import numpy.testing as nt
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


def test_iterable_to_ndarray_for_testing():
    input = [[1, 2, 3], [2, 3, 4]]
    desired = np.array(input)
    actual = ai.iterable_2D_to_ndarray(input)
    nt.assert_allclose(actual=actual, desired=desired, rtol=1e-12)
