
import glob
import re
import sys
import os

import numpy as np

import rdkit
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem

import rmsd

from chemhelp import cheminfo
from chemhelp import misc

# TODO get settings from ini file

global __rungms__
global __tmp__
global __scr__


__rungms__ = "/home/charnley/opt/gamess/gamess-github/rungms"
__scr__ = "/home/charnley/scr/"
__tmp__ = "/home/charnley/scr/"


re_coordinates = re.compile('COORDINATES OF ALL ATOMS ARE (.*?)------------\n(.*?)\n\n', re.DOTALL)
re_error = re.compile('^( \*\*\*|Error:)')

def calculate_optimize(molobj,
    header=None,
    parser=None,
    **kwargs):
    """
    Optimize, get coordinates


    """

    if parser is None:
        parser = read_properties_coordinates

    if header is None:
        header = """ $basis gbasis=pm3 $end
 $contrl runtyp=optimize icharg={:} $end
 $statpt opttol=0.0005 nstep=300 projct=.F. $end
"""

    inpstr = molobj_to_gmsinp(molobj, header)

    stdout, status = run(inpstr, **kwargs)
    properties = parser(stdout)

    return properties


def calculate_thermodynamics(molobj):

    header = """

"""


    return dict()


def calculate_orbitals(molobj):


    return dict()


def calculate_vibrations():


    return dict()


def calculate_solvation():


    return dict()


def calculate(molobj, header, **kwargs):

    inpstr = molobj_to_gmsinp(molobj, header)

    properties = run(inpstr, **kwargs)

    return properties


def prepare_atoms(atoms, coordinates):

    lines = []
    line = "{:2s}    {:2.1f}    {:f}     {:f}    {:f}"

    for atom, coord in zip(atoms, coordinates):
        idx = cheminfo.atom_int(atom)
        lines.append(line.format(atom, idx, *coord))

    lines = [" $data", "Title", "C1"] + lines + [" $end"]

    return "\n".join(lines)



def prepare_xyz(filename, charge, header):
    """
    """

    atoms, coordinates = rmsd.get_coordinates_xyz("test.xyz")

    lines = prepare_atoms(atoms, coordinates)
    header = header.format(charge)

    gmsin = header + lines

    return gmsin



def prepare_mol(filename, header, add_hydrogens=True):
    """
    """

    atoms = []
    coordinates = []

    with open(filename, 'r') as f:
        molfmt = f.read()
        mol = Chem.MolFromMolBlock(molfmt)

    # get formal charge
    charge = Chem.GetFormalCharge(mol)

    # Add hydrogens
    if add_hydrogens:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Get coordinates
    conf = mol.GetConformer(0)
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz = [pos.x, pos.y, pos.z]
        coordinates.append(xyz)
        atoms.append(atom.GetSymbol())

    # set charge
    header = header.format(charge)
    lines = prepare_atoms(atoms, coordinates)

    return header + lines


def molobj_to_gmsinp(mol, header, add_hydrogens=False):
    """
    RDKit Mol object to GAMESS input file

    returns:
        str - GAMESS input file
    """

    coordinates = []
    atoms = []

    # get formal charge
    charge = Chem.GetFormalCharge(mol)

    # Add hydrogens
    if add_hydrogens:
        mol = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol, AllChem.ETKDG())

    # Get coordinates
    conf = mol.GetConformer(0)
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz = [pos.x, pos.y, pos.z]
        coordinates.append(xyz)
        atoms.append(atom.GetSymbol())

    header = header.format(charge)
    lines = prepare_atoms(atoms, coordinates)

    return header + lines


# def calculate(filename, folder="./", store_output=True):
#     """
#     Use GAMESS shell and calculate
#     """
#
#     logfile = filename.replace(".inp", ".log")
#
#     cmd = __rungms__ + " " + filename
#     if store_output:
#         cmd += " > " + folder + logfile
#
#     stdout, stderr = shell.shell(cmd, shell=True)
#
#     return stdout, stderr


def run(inpstr, scr=__tmp__,
    filename=None,
    autoclean=True,
    debug=False):
    """
    """

    if filename is None:
        pid = os.getpid()
        pid = str(pid)
        filename = "_tmp_gamess_run_" + pid + ".inp"

    pwd = os.getcwd()
    os.chdir(scr)

    with open(filename, 'w') as f:
        f.write(inpstr)

    cmd = __rungms__ + " " + filename
    stdout, stderr = misc.shell(cmd)

    if debug:
        print(stderr)

    if autoclean:
        files = glob.glob(__scr__ + filename.replace(".inp", "") + "*")
        for f in files:
            os.remove(f)

    os.chdir(pwd)

    return stdout, stderr


def check_output(output):

    # TODO ELECTRONS, WITH CHARGE ICHARG=

    # TODO redo in Python. Real categories of fail. Better output
    # TODO Did gasphase or solvent phase not converge?
    # TODO How many steps?
    #
    # grep  "FAILURE" *.log
    # grep  "Failed" *.log
    # grep  "ABNORMALLY" *.log
    # grep  "SCF IS UNCONVERGED" *.log
    # grep "resubmit" *.log
    # grep "IMAGINARY FREQUENCY VIBRATION" *.log


    return True, ""


def check_filename(filename):
    """
    Check that the GAMESS calculation is not crashed
    """

    return True, ""


def read_properties(output):

    return

def read_properties_coordinates(output):

    properties = {}

    lines = output.split("\n")

    idx = misc.get_index(lines, "TOTAL NUMBER OF ATOMS")
    line = lines[idx]
    line = line.split("=")
    n_atoms = int(line[-1])

    idx = misc.get_rev_index(lines, "EQUILIBRIUM GEOMETRY LOCATED")
    idx += 4

    coordinates = np.zeros((n_atoms, 3))
    atoms = np.zeros(n_atoms, dtype=int)

    for i in range(n_atoms):
        line = lines[idx + i]
        line = line.split()
        atom = line[1].replace(".0", "")
        atom = int(atom)
        x = line[2]
        y = line[3]
        z = line[4]

        atoms[i] = atom
        coordinates[i][0] = x
        coordinates[i][1] = y
        coordinates[i][2] = z

    idx = misc.get_rev_index(lines, "HEAT OF FORMATION IS")
    line = lines[idx]
    line = line.split()
    hof = float(line[4]) # kcal/mol

    properties["atoms"] = atoms
    properties["coord"] = coordinates
    properties["h"] = hof

    return properties


def read_properties_vibration(output):

    properties = {}

    lines = output.split("\n")

    # Get number of atoms
    idx = misc.get_index(lines, "TOTAL NUMBER OF ATOMS")
    line = lines[idx]
    line = line.split("=")
    n_atoms = int(line[-1])

    # Get heat of formation
    idx = misc.get_rev_index(lines, "HEAT OF FORMATION IS")
    line = lines[idx]
    line = line.split()
    hof = float(line[4]) # kcal/mol

    # Check linear
    idx = misc.get_index(lines, "THIS MOLECULE IS RECOGNIZED AS BEING LINEAR")
    is_linear = (idx is not None)

    # thermodynamic
    idx = misc.get_rev_index(lines, "KJ/MOL    KJ/MOL    KJ/MOL   J/MOL-K")
    idx += 1
    values = np.zeros((5,6))
    for i in range(5):
        line = lines[idx +i]
        line = line.split()
        line = line[1:]
        line = [float(x) for x in line]
        values[i,:] = line

    # Get Vibrations
    idx_start = misc.get_rev_index(lines, "FREQ(CM**-1)")
    idx_end = misc.get_rev_index(lines, "THERMOCHEMISTRY AT T=  298.15 K")
    idx_start += 1
    idx_end -= 2
    vibrations = []
    intensities = []
    for i in range(idx_start, idx_end):
        line = lines[i]
        line = line.split()
        freq = line[1]
        freq = float(line[1])
        inte = line[-1]
        inte = float(inte)
        vibrations.append(freq)
        intensities.append(inte)

    # Cut and save vibration string for jsmol
    # based on number of vibrations and number of atoms
    idx = misc.get_rev_index(lines, " TAKEN AS ROTATIONS AND TRANSLATIONS.")
    vib_lines = "\n".join(lines[idx:idx_start])

    idx_end = misc.get_index(lines, "ELECTRON INTEGRALS")
    head_lines = "\n".join(lines[18:idx_end])

    properties["jsmol"] = head_lines + vib_lines
    properties["linear"] = is_linear
    properties["freq"] = np.array(vibrations)
    properties["intens"] = np.array(intensities)
    properties["thermo"] = values
    properties["h"] = hof

    return properties


def read_properties_orbitals(output):

    properties = {}

    lines = output.split("\n")
    n_lines = len(lines)

    # Get number of atoms
    idx = misc.get_index(lines, "TOTAL NUMBER OF ATOMS")
    line = lines[idx]
    line = line.split("=")
    n_atoms = int(line[-1])

    idx_start = misc.get_index(lines, "EIGENVECTORS")
    idx_start += 4
    idx_end = misc.get_index(lines, "END OF RHF CALCULATION", offset=idx_start)
    energies = []

    wait = False
    j = idx_start
    while j < idx_end:

        line = lines[j].strip()

        if wait:
            if line == "":
                j += 1
                wait = False

        else:
            wait = True

            line = line.split()
            line = [float(x) for x in line]
            energies += line

        j += 1

    properties["orbitals"] = np.array(energies)

    return properties


def read_properties_solvation(output):

    cal_to_joule = 4.18; # cal/mol to joule

    properties = {}

    lines = output.split("\n")
    n_lines = len(lines)

    # Get number of atoms
    idx = misc.get_index(lines, "TOTAL NUMBER OF ATOMS")
    line = lines[idx]
    line = line.split("=")
    n_atoms = int(line[-1])

    # Get solvation data,charge of molecule, surface area, dipole

    idx = misc.get_rev_index(lines, "ELECTROSTATIC INTERACTION")
    line = lines[idx]
    line = line.split()
    electrostatic_interaction = float(line[-2])

    line = lines[idx+1].split()
    pierotti_cavitation_energy = float(line[-2])

    line = lines[idx+2].split()
    dispersion_free_energy = float(line[-2])

    line = lines[idx+3].split()
    repulsion_free_energy = float(line[-2])

    line = lines[idx+4].split()
    total_interaction = float(line[-2])

    total_non_polar = pierotti_cavitation_energy + dispersion_free_energy + repulsion_free_energy


    idx = misc.get_index(lines, "CHARGE OF MOLECULE")
    line = lines[idx]
    line = line.split("=")
    charge = int(line[-1])

    idx = misc.get_rev_index(lines, "SURFACE AREA")
    line = lines[idx]
    line = line.split()
    surface_area = line[2]
    surface_area = surface_area.replace("(A**2)", "")
    surface_area = float(surface_area)

    idx = misc.get_rev_index(lines, "DEBYE")
    line = lines[idx+1]
    line = line.split()
    line = [float(x) for x in line]
    dxyz = line[0:3]
    dtot = line[-1]

    idx = misc.get_rev_index(lines, "MOPAC CHARGES")
    idx += 3
    partial_charges = np.zeros(n_atoms)
    for i in range(n_atoms):
        line = lines[idx+i]
        line = line.split()
        atom_charge = float(line[-2])
        partial_charges[i] = atom_charge

    properties["charges"] = partial_charges
    properties["solvation_total"] = total_interaction
    properties["solvation_polar"] = electrostatic_interaction
    properties["solvation_nonpolar"] = total_non_polar
    properties["surface"] = surface_area
    properties["total_charge"] = charge
    properties["dipole"] = np.array(dxyz)
    properties["dipole_total"] = dtot

    return properties


if __name__ == "__main__":


    with open("header_pm3_opt") as f:
        header = f.readlines()
        header = "".join(header)


    # gmsinp = prepare_xyz("test.xyz", 0, header)

    gmsinp = prepare_mol("test_charge.mol", header)

    f = open("test.inp", 'w')
    f.write(gmsinp)
    f.close()

    calculate("test.inp", folder="./")
    clean() # In production this should be automatic!




