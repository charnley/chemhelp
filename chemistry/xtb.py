
import sys
import os

try:
    import cheminfo

except ModuleNotFoundError:
    here = os.path.dirname(os.path.abspath(__file__))
    sys.path.append(here + "/..")
    import cheminfo

import json
import subprocess
import numpy as np
import ase
import ase.optimize as optlib
import misc

import xtbcalc


def main():


    return


def test():

    atoms = ["C",  "C",  "C",  "C",  "O",  "C",  "O",  "C",  "C",  "H",  "H",  "H",  "H",  "H",  "H",  "H",  "H",  "H",  "H"]

    coord = \
    [[0.449119,  -0.282661,   -0.091832],
    [1.448492,  -0.153837,   0.810153],
    [-0.46586,  -1.171618,   -0.31997],
    [-0.07529,  1.367568,  -2.435109],
    [1.32193,   2.255014,   -3.960735],
    [0.803532,   2.693787,   -2.532549],
    [0.106147,   2.176489,   -0.978965],
    [-0.938151,   1.682795,   -2.615864],
    [-0.305765,   -0.112672,   -1.205421],
    [1.039162,   0.361551,   -0.750651],
    [0.108126,   0.417708,   -0.260873],
    [1.108669,   -2.406275,   0.88672],
    [1.110807,   -1.471904,   0.014561],
    [0.397172,  -2.335684,   -0.651724],
    [1.37253,  1.265402,   -3.871371],
    [0.553584,  3.945833,   -1.63672],
    [-0.042006,  3.5749,   -3.153795],
    [-0.984751,  1.844083,   -0.845294],
    [-1.290445,  1.607219,   -1.111061]]
    coord = np.array(coord)

    atoms_water = ["O", "H", "H"]
    coord_water = [
        [   0.00000  ,  0.00000  ,  0.00000],
        [   0.95094  ,  0.00000  ,  0.00000],
        [   -0.28902 ,  -0.90595 ,  0.00000],
    ]




    molecule = ase.Atoms(atoms, coord)

    calc = xtbcalc.GFN0()

    molecule.set_calculator(calc)
    e = molecule.get_potential_energy()
    print('Initial energy: eV, Eh',e, e/ase.units.Hartree)

    # dyn = LBFGS(molecule)
    # dyn = optlib.LBFGS(molecule)
    dyn = optlib.QuasiNewton(molecule)
    dyn.run(fmax=0.5)

    e = molecule.get_potential_energy()
    print('Initial energy: eV, Eh',e, e/ase.units.Hartree)

    # properties = calculate(atoms_water, coord_water)
    # print(properties)

    return


if __name__ == '__main__':
    test()
