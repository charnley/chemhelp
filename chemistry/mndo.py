
import sys
import os

here = os.path.dirname(__file__)
sys.path.append(here + "/..")

import json
import subprocess
import numpy as np
import ase
import cheminfo
import misc

DEFAULT_PARAMETERS = {
    "method": "PM3",
    "keywords": "precise"
}

MNDOCMD = "mndo.exe"
MNDOCMD = "/home/charnley/opt/mndo/mndo/mndo99_20121112_intel64_ifort-11.1.080_mkl-10.3.12.361"

COORDLINE = "{:} {:} 1 {:} 1 {:} 1"


def get_index(lines, pattern):
    for i, line in enumerate(lines):
        if line.find(pattern) != -1:
            return i
    return None


def reverse_enum(L):
    for index in reversed(range(len(L))):
        yield index, L[index]


def get_rev_index(lines, pattern):

    for i, line in reverse_enum(lines):
        if line.find(pattern) != -1:
            return i

    return None


def run_mndo(stdin):
    """

    just a bunch of notes

    """

    # TODO diskless mndo
    # TODO From procs write "heredocs" << EOF EOF to mndo?

    cmd = MNDOCMD

    stdin = stdin.encode()

    cmd += " << EOF\nPM3 \nEOF"

    print(cmd)
    print()


    # a = "A String of Text"
    # p = subprocess.Popen("wc <<DATA\n" + a + "\nDATA", shell=True)
    # stdout, stderr = p.communicate()


    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    # proc.stdin.close()
    # stdout, stderr = proc.communicate(input="PM3")
    # stdout, stderr = proc.communicate(input=stdin)
    stdout, stderr = proc.communicate()
    # proc.stdin.close()

    print(stdout)

    return stdout, stderr


def run_mndo_file(filename):

    cmd = MNDOCMD
    cmd += " < " + filename

    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stdin=subprocess.PIPE, stderr=subprocess.STDOUT, shell=True)
    stdout, stderr = proc.communicate()

    return stdout, stderr


def create_input(atoms, coordinates, parameters=DEFAULT_PARAMETERS):

    header = []
    header += [parameters["method"]]
    header += [parameters["keywords"]]
    header = " ".join(header)

    lines = []

    for atom, coord in zip(atoms, coordinates):
        lines += [COORDLINE.format(atom, *coord)]

    lines = "\n".join(lines)

    mndoinp = header
    mndoinp += "\n\n\n"
    mndoinp += lines

    return mndoinp


def get_properties(lines):
    """

    """

    # TODO UNABLE TO ACHIEVE SCF CONVERGENCE

    properties = {}

    # Dipole
    idx = get_rev_index(lines, "DIPOLE")
    idx += 3
    line = lines[idx]
    line = line.split()
    value = line[-1]
    value = float(value)
    properties["mu"] = value

    # Enthalpy of formation
    idx_hof = get_index(lines, "FINAL HEAT OF FORMATION")
    line = lines[idx_hof]
    line = line.split("FORMATION")
    line = line[1]
    line = line.split()
    value = line[0]
    value = float(value)
    properties["h"] = value # kcal/mol


    # coordinates
    i = get_rev_index(lines, 'CARTESIAN COORDINATES')
    idx_atm = 1
    idx_x = 2
    idx_y = 3
    idx_z = 4
    n_skip = 4

    if i < idx_hof:
        i = get_rev_index(lines, 'X-COORDINATE')
        idx_atm = 1
        idx_x = 2
        idx_y = 4
        idx_z = 6
        n_skip = 3

    j = i + n_skip
    symbols = []
    coord = []

    # continue until we hit a blank line
    while not lines[j].isspace() and lines[j].strip():
        l = lines[j].split()
        symbols.append(int(l[idx_atm]))
        x = l[idx_x]
        y = l[idx_y]
        z = l[idx_z]
        xyz = [x, y, z]
        xyz = [float(c) for c in xyz]
        coord.append(xyz)
        j += 1

    coord = np.array(coord)
    properties["coord"] = coord
    properties["atoms"] = symbols

    return properties


def read_properties(filename):
    """
    read properties from out file
    """

    with open(filename, 'r') as f:
        lines = f.readlines()

    properties = get_properties(lines)

    return properties


def calculate(atoms, coord, label=None, **kwargs):

    if label is None:
        label = "_tmp_mndo_"

    inpstr = create_input(atoms, coord, **kwargs)

    f = open(label + ".mop", 'w')
    f.write(inpstr)
    f.close()

    stdout, stderr = run_mndo_file(label + ".mop")
    stdout = stdout.decode()
    stdout = stdout.split("\n")

    properties = get_properties(stdout)

    return properties


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

    # inpstr = create_input(atoms, coord)
    # filename = "_tmp_testmndo"
    #
    # f = open(filename + ".mop", 'w')
    # f.write(inpstr)
    # f.close()
    #
    # stdout, stderr = run_mndo_file(filename + ".mop")
    # stdout = stdout.decode()
    #
    # f = open(filename + ".out", 'w')
    # f.write(stdout)
    # f.close()
    #
    # properties = read_properties(filename + ".out")
    # print(properties)

    properties = calculate(atoms_water, coord_water)
    print(properties)

    return


if __name__ == '__main__':
    test()
