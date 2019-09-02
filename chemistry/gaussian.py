
# Add parent

import select
import sys
import os

here = os.path.dirname(__file__)
sys.path.append(here + "/..")

import json
from ase.calculators.gaussian import Gaussian as GaussianCalculator
import cheminfo
import misc
import numpy as np
import ase

import logging
import cclib
from cclib.parser import Gaussian as GaussianReader

import subprocess


DEFAULT_PARAMETERS = {
    "method": "B3LYP",
    "basis": "6-31+",
    "opt": None, # "()"
    "force": None,
    "multiplicity": None
}

def run_gaussian(command):

    errorcode = subprocess.call(command, shell=True)
    print("finished process", errorcode)

    return errorcode

def calculate(atoms, coordinates, parameters=DEFAULT_PARAMETERS, label=None, write_only=True, n_threads=1, mem=1):
    """

    """

    if label is None:
        filename = "_tmp_" + "gaussian"
    else:
        filename = label

    comfile = filename+".com"

    method = parameters["method"]
    basis = parameters["basis"]

    if parameters["opt"] is not None:
       calculator = GaussianCalculator(method=method, basis=basis, opt=parameters["opt"], Freq="()", label=filename)
    else:
       calculator = GaussianCalculator(method=method, basis=basis, label=filename)

    if parameters['force'] is None:
        del calculator.parameters['force']

    molecule = ase.Atoms(atoms, coordinates)
    molecule.set_calculator(calculator)


    calculator.write_input(molecule)

    if n_threads > 1:
        jstr = r"%nprocshared={:}\n%mem={:}GB".format(n_threads, mem)
        cmd = ["sed", "-i", "' 1 s/.*/" + jstr + r"\n&/'", comfile]
        cmd = " ".join(cmd)
        subprocess.call(cmd, shell=True)


    if write_only:
        return True


    # Calculate
    command = calculator.command
    command = command.replace("PREFIX", filename)
    run_gaussian(command)

    # calculator = molecule.get_calculator()
    # calculator.calculate()

    # molecule.get_potential_energy()

    help(calculator)

    # TODO Get All Properties
    properties = read_properties(filename+".log")

    return properties


def read_properties_g4mp2(filename, values={}):


    return


def read_line(filename, pattern):

    for i, line in enumerate(readlines_reverse(filename)):
        if line.find(pattern) != -1:
            return line

    return None


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


def read_properties_b3lyp(filename, values={}):
    """
    """

    # gausslog = GaussianReader(filename)
    # gausslog.logger.setLevel(logging.ERROR)
    # data = gausslog.parse()

    # help(data)
    # data.writejson(filename="_tmp_.json")


    properties = {}

    #  3  A         GHz          Rotational constant A
    #  4  B         GHz          Rotational constant B
    #  5  C         GHz          Rotational constant C
    #  6  mu        Debye        Dipole moment
    #  7  alpha     Bohr^3       Isotropic polarizability
    #  8  homo      Hartree      Energy of Highest occupied molecular orbital (HOMO)
    #  9  lumo      Hartree      Energy of Lowest occupied molecular orbital (LUMO)
    # 10  gap       Hartree      Gap, difference between LUMO and HOMO
    # 11  r2        Bohr^2       Electronic spatial extent
    # 12  zpve      Hartree      Zero point vibrational energy
    # 13  U0        Hartree      Internal energy at 0 K
    # 14  U         Hartree      Internal energy at 298.15 K
    # 15  H         Hartree      Enthalpy at 298.15 K
    # 16  G         Hartree      Free energy at 298.15 K
    # 17  Cv        cal/(mol K)  Heat capacity at 298.15 K

    f = open(filename, 'r')

    lines = f.readlines()

    # Rotational constants
    idx = get_rev_index(lines, "Rotational constants")
    rots = lines[idx].split()
    properties["A"] = float(rots[-3])
    properties["B"] = float(rots[-2])
    properties["C"] = float(rots[-1])

    # Dipole moment
    #wrong idx = get_rev_index(lines, "Dipole        =")
    #wrong line = lines[idx].replace("D", "E")
    #wrong offset = 16
    #wrong length = 15
    #wrong value1 = line[offset:offset+length]
    #wrong value2 = line[offset+length:offset+length*2]
    #wrong value3 = line[offset+length*2:offset+length*3]
    idx = get_rev_index(lines, "Dipole moment (field-independent basis, Debye)")
    line = lines[idx+1].split()
    # dipole = [line[1], line[3], line[5]]
    total_dipole = float(line[-1])
    properties["mu"] = total_dipole

    # Isotropic polarizability
    idx = get_rev_index(lines, "Isotropic polarizability")
    line = lines[idx]
    line = line.split()
    properties["alpha"] = float(line[-2])

    # homo and lumo
    idx = get_rev_index(lines, "The electronic state is")
    homo, lumo = get_moenergies(lines, idx+1)
    properties["homo"] = homo
    properties["lumo"] = lumo
    properties["gap"] = lumo - homo

    # Electronic spatial extent
    idx = get_rev_index(lines, "Electronic spatial extent")
    line = lines[idx]
    line = line.split()
    properties["r2"] = float(line[-1])

    # zpve
    idx = get_rev_index(lines, "Zero-point correction")
    line = lines[idx]
    line = line.split()
    properties["zpve"] = float(line[-2])

    # u0
    idx = get_rev_index(lines, "Sum of electronic and zero-point Energies")
    line = lines[idx]
    line = line.split()
    properties["u0"] = float(line[-1])

    # U
    idx = get_rev_index(lines, "Sum of electronic and thermal Energies")
    line = lines[idx]
    line = line.split()
    properties["u"] = float(line[-1])

    # H
    idx = get_rev_index(lines, "Sum of electronic and thermal Enthalpies")
    line = lines[idx]
    line = line.split()
    properties["h"] = float(line[-1])

    # G
    idx = get_rev_index(lines, "Sum of electronic and thermal Free Energies")
    line = lines[idx]
    line = line.split()
    properties["g"] = float(line[-1])

    # Cv
    idx = get_rev_index(lines, "Cal/Mol-Kelvin")
    idx += 1
    line = lines[idx]
    line = line.split()
    cv = line[2]
    properties["cv"] = float(cv)

    return properties


def get_moenergies(lines, idx_start):

    homos = []
    lumos = []

    i = idx_start
    line = lines[i]

    while line.find('Alpha') == 1:

        values = line.split("--")

        energies = values[-1]
        energies = energies.split()
        energies = [float(x) for x in energies]

        if line.find("virt") != -1:
            lumos += energies
        else:
            homos += energies

        # next
        i += 1
        line = lines[i]

    homo = max(homos)
    lumo = min(lumos)

    return homo, lumo



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


def calculate_energy_atomization(atoms, coordinates, method="B3LYP", basis="6-31+", atom_energies=None):

    if atom_energies is None:
        atom_energies = load_atom_energies(basis=basis)

    calculator = Gaussian(method=method, basis=basis, opt="()")
    del calculator.parameters['force']

    e_atom = 0
    for atom in atoms:
        energy = atom_energies[atom]
        e_atom += energy

    molecule = ase.Atoms(atoms, coordinates)
    molecule.set_calculator(calculator)

    e_molecule = molecule.get_potential_energy()
    e_atomization = e_molecule - 2 * e_atom

    return -e_atomization


def get_atomization(method, basis):

    energies = {}
    atoms = cheminfo.MULTIPLICITY.keys()

    for atom in atoms:

        multiplicity = cheminfo.MULTIPLICITY[atom]

        label = "_tmp_atom_"+atom

        calculator = GaussianCalculator(
            method=method,
            basis=basis,
            multiplicity=multiplicity,
            label=label)

        del calculator.parameters['force']

        atomobj = ase.Atoms(atom, calculator=calculator)

        calculator.write_input(atomobj)
        command = calculator.command
        command = command.replace("PREFIX", label)
        run_gaussian(command)

        if "G4" in method:
            pattern = method + " Enthalpy"
            line = read_line(label + ".log", pattern)
            line = line.split()
            value = line[2]
        else:
            line = read_line(label + ".log", "SCF Done:")
            line = line.split()
            value = line[4]

        value = float(value)
        energies[atom] = value

        print(atom, value)

    return energies


def atomization(atoms, atom_energies, e_molecule):

    e_atom = 0
    for atom in atoms:
        energy = atom_energies[atom]
        e_atom += energy

    e_atomization = e_molecule - e_atom

    return e_atomization


def calculate_energy_atomization(atoms, coordinates, method="B3LYP", basis="6-31+", atom_energies=None):

    if atom_energies is None:
        atom_energies = load_atom_energies(basis=basis)

    calculator = GaussianCalculator(method=method, basis=basis, opt="()")
    del calculator.parameters['force']

    e_atom = 0
    for atom in atoms:
        energy = atom_energies[atom]
        e_atom += energy

    molecule = ase.Atoms(atoms, coordinates)
    molecule.set_calculator(calculator)

    e_molecule = molecule.get_potential_energy()
    e_atomization = e_molecule - 2 * e_atom

    return -e_atomization



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
        name = name.replace(".log", "")

        try:
            properties = read_properties_b3lyp(line)
            if debug: print(name)
            yield name, properties
        except:
            misc.eprint(name, "fail")
            continue


def readlines_sdf():
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
        name = name.replace(".sdf", "")

        with open(line, 'r') as f:
            sdf = f.read()
            yield name, sdf


def print_csv(data, sep=", "):

    data = [str(x) for x in data]

    line = sep.join(data)

    return line


def qm9_format(name, properties):

    data = []

    keys = ['A', 'B', 'C', 'mu', 'alpha', 'homo', 'lumo', 'gap', 'r2', 'zpve', 'U0', 'U', 'H', 'G', 'cv']


    for key in keys:
        data.append(properties[key])

    data = [name] + data
    data = [str(x) for x in data]
    data = ", ".join(data)

    return data


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('--log', action='store', help='', metavar='FILE')
    parser.add_argument('--stdlog', action='store_true', help='')
    parser.add_argument('--stdsdf', action='store_true', help='')
    parser.add_argument('--atom', action='store_true', help='atomization correction')

    parser.add_argument('--dryrun', action='store_true', help='only write input files')
    parser.add_argument('--prefix', action='store', help='', metavar='STR', default="")

    parser.add_argument('-p', '--parameters', help='Parameters for gaussian calculation')
    parser.add_argument('-j', '--threads', help="Prepare multi-thread gaussian jobs", metavar="int", default=0)

    args = parser.parse_args()


    if args.parameters is None:
        parameters = {
            "method": "B3LYP",
            "basis": "6-31G(2df,p)",
            "opt": "()",
            "force": None }
        s = json.dumps(parameters)
        print(s)

    else:
        with open(args.parameters, 'r') as f:
            dumpstr = f.read()
        parameters = json.loads(dumpstr)


    # if args.stdsdf:
    #     for name, sdf in readlines_sdf():
    #
    #         print(name)
    #
    #         molobj, status = cheminfo.sdfstr_to_molobj(sdf)
    #
    #         atoms, coordinates = cheminfo.molobj_to_xyz(molobj)
    #         atoms_str = [cheminfo.convert_atom(atom) for atom in atoms]
    #
    #         properties = calculate(atoms, coordinates, parameters=parameters, label="jobs/"+name, n_threads=6, mem=14)


    if args.sdf:

        with open(args.sdf, 'r') as f:
            sdf = f.read()

        filename = args.sdf
        filename = filename.replace(".sdf", "")

        molobj, status = cheminfo.sdfstr_to_molobj(sdf)

        atoms, coordinates = cheminfo.molobj_to_xyz(molobj)
        atoms_str = [cheminfo.convert_atom(atom) for atom in atoms]

        methodname = parameters["method"]
        basis = parameters["basis"]

        if args.atom:
            energies = get_atomization(methodname, basis)
            corr = atomization(atoms_str, energies, 0.0)
            print("correction", corr)


        else:
            properties = calculate(atoms, coordinates, parameters=parameters, write_only=True, label=filename+args.prefix, n_threads=30, mem=150)
            pass

        # e_atomi = calculate_energy_atomization(atoms_str, coordinates, atom_energies=energies, method="B3LYP", basis="6-31G(2df,p)")
        #
        # print(e_atomi)


    if args.stdlog:

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




    if args.log:

        properties = read_properties(args.log)

        for key in properties.keys():

            fmt = "{:8s}          {:15.8f}".format(key, properties[key])
            print(fmt)

    return

if __name__ == '__main__':
    main()

