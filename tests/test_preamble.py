import numpy as np


def atom_equal(at1, at2):
    if len(at1) != len(at2):  # pragma: no branch
        return False
    elif at1[0] != at2[0]:
        return False
    elif not np.allclose(at1[1], at2[1], atol=1.e-4):
        return False
    elif at1[2:] != at2[2:]:
        return False
    else:
        return True


def vec_equal(v1, v2):
    return np.allclose(v1, v2, atol=1.e-6)


def float_equal(f1, f2):
    return np.isclose(f1, f2, atol=1.e-6)
