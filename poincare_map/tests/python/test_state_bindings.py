import numpy as np
import numpy.testing as nt
import dynamic_analysis.core as dyn_core


def test_args_to_array():
    result = dyn_core.args_to_array(1., 2., 3.)
    desired = np.array([1., 2., 3.])
    nt.assert_array_equal(result, desired)


def test_args_to_Hessian():
    result = dyn_core.args_to_Hessian(1, 2, 3)
    desired = np.array([[1., 2.],
                        [2., 3.]])
    nt.assert_array_equal(result, desired)


def test_arma_mat_copying():
    input = [[1, 2, 3], [2, 3, 4]]
    desired = np.array(input)
    actual = dyn_core.iterable_2D_to_ndarray(input)
    nt.assert_allclose(actual=actual, desired=desired, rtol=1e-12)


def test_std_vector_double_copying():
    input = [1, 2, 3, 2, 3, 4]
    desired = np.array(input)
    actual = dyn_core.iterable_1D_to_ndarray(input)
    nt.assert_allclose(actual=actual, desired=desired, rtol=1e-12)
