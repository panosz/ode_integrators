import pytest
from numpy.testing import assert_array_equal
import dynamic_analysis.near_resonance as nr


def test_Resonance_invalid_string_construction():
    invalid_string = '1 over 3'
    with pytest.raises(ValueError) as constrinfo:
        nr.Resonance(invalid_string)
    assert 'Invalid literal' in str(constrinfo.value)


def test_Resonance_invalid_type_construction():
    invalid_input = (1, 2, 3)
    with pytest.raises(TypeError) as constrinfo:
        nr.Resonance(invalid_input)
    assert 'Invalid resonance' in str(constrinfo.value)


def test_Resonance_valid_string_construction():
    res = nr.Resonance("1/2")
    assert res.ratio == 0.5


def test_Resonance_valid_iterable_construction():
    res = nr.Resonance([1, 2])
    assert res.ratio == 0.5


def test_Resonance_zero_denominator_exception():
    with pytest.raises(ZeroDivisionError):
        nr.Resonance((1, 0))


def test_Resonance_ratio_read_only():
    res = nr.Resonance([1, 2])
    with pytest.raises(AttributeError) as atrerror:
        res.ratio = 3
    assert "can't set attribute" in str(atrerror.value)


def test_Resonance_jacobian():
    num = 3
    den = 4
    res = nr.Resonance((num, den))
    desired = [[num, 0],
               [-den, 1]]
    assert_array_equal(res.jacobian(), desired)


def test_Resonance_transform_hessian():
    res = nr.Resonance('3/5')
    input_hess = [[2, 4], [4, -7]]
    desired = [[-277, 47], [47, -7]]
    actual = res.transform_hessian(input_hess)
    assert_array_equal(actual, desired)


def test_Resonance_comparison():
    assert nr.Resonance((3, 6)) == nr.Resonance((-1, -2))
    assert nr.Resonance(3, 6) == 1/2


def test_Resonance_repr():
    assert "Resonance(1, 2)" in repr(nr.Resonance(2, 4))


def test_Resonance_str():
    assert "1/2" in str(nr.Resonance(2, 4))


def test_not_addable():
    with pytest.raises(TypeError):
        nr.Resonance(2, 4) + nr.Resonance(3, 5)


def test_comparable():
    assert nr.Resonance(3, 6) > nr.Resonance(3, 7)
    assert nr.Resonance(3, 6) >= 0.5
    assert nr.Resonance(3, 6) <= 0.5
    assert nr.Resonance(-3, 6) < 0
