# -*- coding: utf-8 -*-

import numpy as np
from fractions import gcd

from vipster import Molecule

Agraph = 4.6631  # bohr
sq3 = np.sqrt(3)
vecX = np.array([1, 0, 0], 'f')
vecY = np.array([0, 1, 0], 'f')
vecZ = np.array([0, 0, 1], 'f')
vecYh = np.array([-0.5, 0.5 * sq3, 0], 'f')


def grapheneHex(x, y, z):
    """
    Create a hexagonal cell of graphene

    x -> x-dimension
    y -> 0.5*x+0.5*sqrt(3)*y-dimension
    z -> interlayer distance
    """
    mol = Molecule(name="({:},{:})-hex graphene".format(x, y))
    mol.setCellDim(Agraph)
    mol.setVec([vecX, vecYh, z * vecZ])
    mol.newAtom('C', [1. / 3., 2. / 3., 0])
    mol.newAtom('C', [2. / 3., 1. / 3., 0])
#    mol.scaleAtoms("crystal")
#    mol.setFmt("crystal")
    mol.setFmt("crystal", scale=True)
    mol.mult(x, y, 1)
    return mol


def grapheneOrtho(x, y, z):
    """
    Create an orthogonal cell of graphene

    x -> x-dimension
    y -> sqrt(3)*y-dimension
    z -> interlayer distance
    """
    mol = Molecule(name="({:},{:})-ortho graphene".format(x, y))
    mol.setCellDim(Agraph)
    mol.setVec([vecX, sq3 * vecY, z * vecZ])
    mol.newAtom('C', [0, 0, 0])
    mol.newAtom('C', [0, 1. / 3., 0])
    mol.newAtom('C', [0.5, 0.5, 0])
    mol.newAtom('C', [0.5, 5. / 6., 0])
#    mol.scaleAtoms("crystal")
#    mol.setFmt("crystal")
    mol.setFmt("crystal", scale=True)
    mol.mult(x, y, 1)
    return mol


def cnt(m, n, l=0):
    """
    Create a carbon nanotube

    m, n -> vector for "roll-up"
    l -> target length in angstrom
    """
    mol = grapheneOrtho(3 * max(m, n), 2 * max(m, n), 1)
    mol.name = "({:},{:}) cnt".format(m, n)
    middle = 0.5 * mol.getCellDim() * mol.getVec()[0]
    mol.shift(range(mol.nat), -middle)
    circumference = np.sqrt(m**2 + n**2 + m * n)
    angle = -np.degrees(np.arccos((2 * m + n) / (2 * circumference)))
    divisor = gcd(2 * m + n, 2 * n + m)
    p1 = (2 * n + m) / divisor
    p2 = -(2 * m + n) / divisor
    periodicity = np.sqrt(p1**2 + p2**2 + p1 * p2)
    mol.rotate(range(mol.nat), angle, [0, 0, 1])
    mol.setVec([circumference * vecX, periodicity * vecY,
                circumference / 2 / np.pi * vecZ])
    mol.crop()
    if True:
        # move above,  then wrap around z-axis
        for i in range(mol.nat):
            at = mol.getAtom(i, fmt="crystal")[1]
            angle = at[0] * 360
            mol.setAtom(i, coord=[0, at[1], 1], fmt="crystal")
            mol.rotate([i], angle, [0, 1, 0])
    # approximate target length
    l = l / (mol.getCellDim(fmt="angstrom") * mol.getVec()[1][1])
    mol.mult(1, max(1, int(l)), 1)
    vec = mol.getVec()
    vec[2][2] = circumference
    mol.setVec(vec)
    return mol
