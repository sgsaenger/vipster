import numpy as np


def atom_equal(at1, at2):
    if len(at1) != len(at2):
        return False
    elif at1[0] != at2[0]:
        return False
    elif not np.allclose(at1[1], at2[1]):
        return False
    elif at1[2:] != at2[2:]:
        return False
    else:
        return True
