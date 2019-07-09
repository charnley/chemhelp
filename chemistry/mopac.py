
import select
import sys
import os

here = os.path.dirname(__file__)
sys.path.append(here + "/..")

import subprocess
import numpy as np
import ase
from ase.calculators.mopac import MOPAC as Mopac

import cheminfo
import misc


def readlines_reverse(filename):
    with open(filename) as qfile:
        qfile.seek(0, os.SEEK_END)
        position = qfile.tell()
        line = ''
        while position >= 0:
            qfile.seek(position)
            next_char = qfile.read(1)
            if next_char == "\n":
                yield line[::-1]
                line = ''
            else:
                line += next_char
            position -= 1
        yield line[::-1]

def read_line(filename, pattern):

    for i, line in enumerate(readlines_reverse(filename)):
        if line.find(pattern) != -1:
            return line

    return None


def run_mopac(command):

    errorcode = subprocess.call(command, shell=True)

    return errorcode


def read_coord_from_mopac_file(filename, find_energy=False):

    with open(filename, 'r') as f:
        lines = f.readlines()

    i = get_index(lines, 'CARTESIAN COORDINATES')

    j = i + 4
    symbols = []
    coord = []

    while not lines[j].isspace():  # continue until we hit a blank line
        l = lines[j].split()
        symbols.append(l[1])
        coord.append([float(c) for c in l[2: 2 + 3]])
        j += 1

    coord = np.array(coord)


    if find_energy:

        i = get_index(lines, 'FINAL HEAT OF FORMATION')

        if i is None:
            energy = None

        else:

            energy = lines[i]
            energy = energy.split("=")[1]
            energy = energy.replace("KCAL/MOL", "")
            energy = float(energy)

        return energy, symbols, coord


    return symbols, coord


def read_properties(filename):

    properties = {}

    with open(filename, 'r') as f:
        lines = f.readlines()

    # Dipole
    idx = get_rev_index(lines, "DIPOLE")
    idx += 3
    line = lines[idx]
    line = line.split()
    value = line[-1]
    value = float(value)
    properties["mu"] = value

    # Enthalpy of formation
    idx = get_index(lines, "FINAL HEAT OF FORMATION")
    line = lines[idx]
    line = line.split("=")
    line = line[1]
    line = line.split()
    value = line[0]
    value = float(value)
    properties["h"] = value

    # homo and lumo
    idx = get_index(lines, "HOMO LUMO ENERGIES")
    line = lines[idx]
    line = line.split()
    homo = float(line[-2]) # ev
    lumo = float(line[-1]) # ev

    homo *= 23.06035 # ev to kcal/mol
    lumo *= 23.06035 # ev to kcal/mol

    properties["homo"] = homo
    properties["lumo"] = lumo
    properties["gap"] = lumo - homo

    return properties


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


def optimize_and_energy_mopac(atoms, coord, filename="", method="PM7"):

    filename += "_opten"

    calculator = Mopac(method=method, task="precise", label=filename)

    molecule = ase.Atoms(atoms, coord)
    molecule.set_calculator(calculator)

    energy = molecule.get_potential_energy()

    return energy


def optimize_mopac(atoms, coord, filename="", method="PM7"):

    filename += "_optimize"

    calculator = Mopac(method=method, task="precise", label=filename)

    molecule = ase.Atoms(atoms, coord)
    molecule.set_calculator(calculator)

    energy = molecule.get_potential_energy()

    symbols, coord = read_coord_from_mopac_file(filename + ".out")

    return coord


def readlines_log(debug=False):
    """
    read stdin, if any.
    Assume sdf file names and return sdfs
    """

    while sys.stdin in select.select([sys.stdin], [], [], 0)[0]:

        line = sys.stdin.readline()

        if not line:
            yield from []
            break

        line = line.strip()

        name = line.split("/")
        name = name[-1]
        name = name.replace(".out", "")

        try:
            properties = read_properties(line)
            if debug: print(name)
            yield name, properties
        except:
            print(name, "fail")
            continue


def get_atomization(method):

    # singlet - 0 unpaired electrons
    # doublet - 1 unpaired electrons
    # triplet - 2 unpaired electrons
    # quartet - 3 unpaired electrons
    # quintet - 4 unpaired electrons
    # sextet  - 5 unpaired electrons

    energies = {}
    atoms = cheminfo.MULTIPLICITY.keys()

    for atom in atoms:

        label = "_tmp_atom_"+atom
        multiplicity = cheminfo.MULTIPLICITY[atom]
        keyword = None

        if multiplicity == 1:
            keyword = "singlet"
        elif multiplicity == 2:
            keyword = "doublet"
        elif multiplicity == 3:
            keyword = "triplet"
        elif multiplicity == 4:
            keyword = "quartet"
        elif multiplicity == 5:
            keyword = "quintet"
        elif multiplicity == 6:
            keyword = "sextet"

        calculator = Mopac(method=method, task="uhf precise charge=0 {}".format(keyword), label=label)
        atomobj = ase.Atoms(atom, calculator=calculator)
        calculator.write_input(atomobj)

        cmd = calculator.command
        cmd = cmd.replace("PREFIX", label)
        run_mopac(cmd)

        # Line
        line = read_line(label + ".out", "FINAL HEAT OF FORMATION")
        print(line, atom)
        line = line.split()

        hof = line[5]
        energies[atom] = float(hof)

    return energies


def atomization(atoms, atom_energies, e_molecule):

    e_atom = 0
    for atom in atoms:
        energy = atom_energies[atom]
        e_atom += energy

    e_atomization = e_molecule - e_atom

    return e_atomization


def print_csv(data, sep=", "):

    data = [str(x) for x in data]

    line = sep.join(data)

    return line


def main():
    import sys
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('--log', action='store', help='', metavar='FILE')
    parser.add_argument('--stdlog', action='store_true', help='')
    parser.add_argument('--stdsdf', action='store_true', help='')
    parser.add_argument('--atom', action='store_true', help='atomization')

    args = parser.parse_args()


    if args.sdf:

        with open(args.sdf, 'r') as f:
            sdf = f.read()

        molobj, status = cheminfo.sdfstr_to_molobj(sdf)

        atoms, coordinates = cheminfo.molobj_to_xyz(molobj)
        atoms_str = [cheminfo.convert_atom(atom) for atom in atoms]

        if args.atom:
            energies = get_atomization("PM6")
            print(atoms_str)
            corr = atomization(atoms_str, energies, 0.0)
            print("correction", corr)


    if args.log:

        properties = read_properties(args.log)

        for key in properties.keys():

            fmt = "{:8s}          {:15.8f}".format(key, properties[key])
            print(fmt)


    if args.stdlog:
        """
        print out qm9 properties

        """

        datout = []

        sep = ", "

        generator = readlines_log(debug=False)

        name, properties = next(generator)

        keys = properties.keys()
        keys = ['mu', 'h', 'homo', 'lumo', 'gap']

        header = print_csv(["name"]+keys)
        print(header)

        data = [name]+[properties[key] for key in keys]
        line = print_csv(data)
        print(line)

        for name, properties in generator:

            data = [name]+[properties[key] for key in keys]
            line = print_csv(data)
            print(line)


    return

if __name__ == '__main__':
    main()
