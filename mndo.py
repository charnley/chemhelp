
from pathlib import Path
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
import subprocess
import json
import os
import functools
import shutil
import copy

from . import units
from . import misc

IGNORE_PARAMS = [
    'DD2',
    'DD3',
    'PO1',
    'PO2',
    'PO3',
    'PO9',
    'HYF',
    'CORE',
    'EISOL',
    'FN1',
    'FN2',
    'FN3',
    'GSCAL',
    'BETAS',
    'ZS',
]

MNDO_CMD = "mndo"

PARAMETERS_OM2 = {
    "H": {
        "USS": -11.90627600,
        "ZS": 1.33196700,
        "BETAS": -6.98906400,
        "ALP": 2.54413410,
    },
    "C": {
        "ZP": 1.78753700,
        "BETAS": -18.98504400,
        "BETAP": -7.93412200,
        "ALP": 2.54638000,
    },
    "N": {
        "USS": -71.93212200,
        "UPP": -57.17231900,
        #"ZS": 2.25561400,
        "ZP": 2.25561400,
        #"BETAS": -20.49575800,
        "BETAP": -20.49575800,
    },
    "O": {
        "USS": -99.64430900,
        "UPP": -77.79747200,
        # "ZS": 2.69990500,
        "ZP": 2.69990500,
        # "BETAS": -32.68808200,
        "BETAP": -32.68808200,
        "ALP": 3.16060400,
    },
}

DEFAULT_HEADER = """{method} 1SCF MULLIK PRECISE charge={charge} jprint=5
nextmol=-1
TITLE {title}"""

ATOMLINE = "{:2s} {:} 0 {:} 0 {:} 0"

__DEBUG__ = "jprint=7"


__ATOMS__ = [x.strip() for x in [
    'h ', 'he', \
    'li', 'be', 'b ', 'c ', 'n ', 'o ', 'f ', 'ne', \
    'na', 'mg', 'al', 'si', 'p ', 's ', 'cl', 'ar', \
    'k ', 'ca', 'sc', 'ti', 'v ', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', \
    'zn', 'ga', 'ge', 'as', 'se', 'br', 'kr',  \
    'rb', 'sr', 'y ', 'zr', 'nb', 'mo', 'tc', 'ru', 'rh', 'pd', 'ag', \
    'cd', 'in', 'sn', 'sb', 'te', 'i ', 'xe',  \
    'cs', 'ba', 'la', 'ce', 'pr', 'nd', 'pm', 'sm', 'eu', 'gd', 'tb', 'dy', \
    'ho', 'er', 'tm', 'yb', 'lu', 'hf', 'ta', 'w ', 're', 'os', 'ir', 'pt', \
    'au', 'hg', 'tl', 'pb', 'bi', 'po', 'at', 'rn', \
    'fr', 'ra', 'ac', 'th', 'pa', 'u ', 'np', 'pu']]


def convert_atom(atom, t=None):
    """

    convert atom from str2int or int2str

    """

    if t is None:
        t = type(atom)
        t = str(t)

    if "str" in t:
        atom = atom.lower()
        idx = __ATOMS__.index(atom) + 1
        return idx

    else:
        atom = __ATOMS__[atom -1].capitalize()
        return atom


def execute(cmd, scr=None):

    if scr is not None:
        cmd = f"cd {scr}; " + cmd

    popen = subprocess.Popen(cmd,
        stdout=subprocess.PIPE,
        universal_newlines=True,
        shell=True)

    for stdout_line in iter(popen.stdout.readline, ""):
        yield stdout_line

    popen.stdout.close()
    return_code = popen.wait()

    if return_code:
        raise subprocess.CalledProcessError(return_code, cmd)


def run_mndo_file(filename, scr=None, mndo_cmd=MNDO_CMD):

    runcmd = f"{mndo_cmd} < {filename}"

    lines = execute(runcmd, scr=scr)

    molecule_lines = []

    for line in lines:

        molecule_lines.append(line.strip("\n"))

        if "STATISTICS FOR RUNS WITH MANY MOLECULES" in line:
            return

        if "COMPUTATION TIME" in line:
            yield molecule_lines
            molecule_lines = []

    return


def calculate_file(filename,
    return_properties=True,
    **kwargs):
    """

    Run MNDO.exe on filename and return iterator of results.

    if return_properties is set, return Iterator[Dict] with properties
    else returns Iterator[List[Str]] of output

    :param filename: 
    :param return_properties: boolean
    :return results: Iterator
    """

    calculations = run_mndo_file(filename, **kwargs)

    if not return_properties:

        if return_list:
            calculations = list(calculations)

        return calculations

    for lines in calculations:
        properties = get_properties(lines)
        yield properties

    return


# def calculate(filename, **kwargs):
#     """
#
#     DEPRECATIED
#
#     TODO rewrite calculate interface
#
#     
#
#     """
#
#     calculations = run_mndo_file(filename, **kwargs)
#
#     properties_list = []
#
#     for lines in calculations:
#         properties = get_properties(lines)
#         properties_list.append(properties)
#
#     return properties_list



def get_properties(lines):
    """

    get properties of a single calculation

    arguments:
        lines - list of MNDO output lines

    return:
        dict of properties

    """

    # TODO Check on input type of calculation

    # TODO UNABLE TO ACHIEVE SCF CONVERGENCE
    # TODO optimization failed

    # TODO Return NaN if idx is not found

    # TODO Check CORE HAMILTONIAN MATR if not converged?

    # check for error
    # idx = get_index(lines, "UNABLE TO ACHIEVE SCF CONVERGENCE")
    # if idx is not None:
    #     properties["energy"] = np.nan
    #     return properties


    properties = get_properties_1scf(lines)

    return properties


def get_properties_1scf(lines):

    properties = {}

    # Check if input coordiantes is internal
    # INPUT IN INTERNAL COORDINATES
    # INPUT IN CARTESIAN COORDINATES
    idx = misc.get_index(lines, "INPUT IN")
    line = lines[idx]
    is_internal = "INTERNAL" in line

    keywords = [
        "CORE HAMILTONIAN MATRIX.",
        "NUCLEAR ENERGY",
        "IONIZATION ENERGY",
        "INPUT GEOMETRY"]

    idx_keywords = misc.get_rev_indexes(lines, keywords)

    # SCF energy
    idx_core = idx_keywords[0]
    if idx_core is None:
        e_scf = float("nan")
        properties["e_scf"] = e_scf
    else:
        idx = idx_core
        idx -= 9
        line = lines[idx]

        if "SCF CONVERGENCE HAS BEE" in line:
            idx -= 2
            line = lines[idx]

        # NOTE This should never happen, but better safe than sorry
        line = line.split()
        if len(line) < 2:
            e_scf = float("nan")
        else:
            value = line[1]
            e_scf = float(value)

        properties["e_scf"] = e_scf # ev

    # Nuclear energy
    if idx_keywords[1] is None:
        e_nuc = float("nan")
        properties["e_nuc"] = e_nuc
    else:
        idx = idx_keywords[1]
        line = lines[idx]
        line = line.split()
        value = line[2]
        e_nuc = float(value)
        properties["e_nuc"] = e_nuc # ev

    # eisol
    eisol = dict()
    idxs = misc.get_indexes_with_stop(lines, "EISOL", "IDENTIFICATION")
    for idx in idxs:
        line = lines[idx]
        line = line.split()
        atom = int(line[0])
        value = line[2]
        eisol[atom] = float(value) # ev

    # # Enthalpy of formation
    idx_hof = misc.get_index(lines, "SCF HEAT OF FORMATION")
    if idx_hof is None:
        properties["h"] = float("nan")
    else:
        line = lines[idx_hof]
        line = line.split("FORMATION")
        line = line[1]
        line = line.split()
        value = line[0]
        value = float(value)
        properties["h"] = value # kcal/mol

    # ionization
    # idx = get_rev_index(lines, "IONIZATION ENERGY")
    idx = idx_keywords[2]
    if idx is None:
        e_ion = float("nan")
        properties["e_ion"] = e_ion
    else:
        line = lines[idx]
        value = line.split()[-2]
        e_ion = float(value) # ev
        properties["e_ion"] = e_ion

    # # Dipole
    # idx = get_rev_index(lines, "PRINCIPAL AXIS")
    # line = lines[idx]
    # line = line.split()
    # value = line[-1]
    # value = float(value) # Debye
    # properties["mu"] = value

    # input coords

    atoms = []
    coord = []

    if is_internal:

        idx_atm = 1
        idx_x = 2
        idx_y = 3
        idx_z = 4

        idx_coord = misc.get_index(lines, "INITIAL CARTESIAN COORDINATES")
        idx_coord += 5

        j = idx_coord
        # continue until we hit a blank line
        while not lines[j].isspace() and lines[j].strip():
            line = lines[j].split()

            atom = line[idx_atm]
            atom = int(atom)
            x = float(line[idx_x])
            y = float(line[idx_y])
            z = float(line[idx_z])

            atoms.append(atom)
            coord.append([x,y,z])

            j += 1

    else:

        idx_atm = 1
        idx_x = 2
        idx_y = 3
        idx_z = 4

        idx = idx_keywords[3]
        idx += 6

        j = idx
        # continue until we hit a blank line
        while not lines[j].isspace() and lines[j].strip():
            l = lines[j].split()
            atoms.append(int(l[idx_atm]))
            x = l[idx_x]
            y = l[idx_y]
            z = l[idx_z]
            xyz = [x, y, z]
            xyz = [float(c) for c in xyz]
            coord.append(xyz)
            j += 1

    # calculate energy
    e_iso = [eisol[a] for a in atoms]
    e_iso = np.sum(e_iso)
    energy = (e_nuc + e_scf - e_iso)

    properties["energy"] = energy

    return properties


def get_properties_optimize():

    # # optimized coordinates
    # i = get_rev_index(lines, 'CARTESIAN COORDINATES')
    # idx_atm = 1
    # idx_x = 2
    # idx_y = 3
    # idx_z = 4
    # n_skip = 4
    #
    # if i < idx_hof:
    #     i = get_rev_index(lines, 'X-COORDINATE')
    #     idx_atm = 1
    #     idx_x = 2
    #     idx_y = 4
    #     idx_z = 6
    #     n_skip = 3
    #
    # j = i + n_skip
    # symbols = []
    # coord = []
    #
    # # continue until we hit a blank line
    # while not lines[j].isspace() and lines[j].strip():
    #     l = lines[j].split()
    #     symbols.append(int(l[idx_atm]))
    #     x = l[idx_x]
    #     y = l[idx_y]
    #     z = l[idx_z]
    #     xyz = [x, y, z]
    #     xyz = [float(c) for c in xyz]
    #     coord.append(xyz)
    #     j += 1
    #
    # coord = np.array(coord)
    # properties["coord"] = coord
    # properties["atoms"] = symbols


    return

def get_properties_gradient():


    return


def _dump_lines(filename, lines):
    """
    """

    if type(lines) == list:
        lines = "\n".join(lines)

    with open(filename, 'w') as f:
        f.write(lines)

    return


def param_worker(*args, **kwargs):

    scr = kwargs["scr"]
    filename = kwargs["filename"]
    params = args[0]

    # Find subfolder
    current = mp.current_process()
    pid = current.name
    pid, = current._identity
    pid = str(pid)
    scr = os.path.join(scr, pid)

    # Set params in worker dir
    set_params(params, scr=scr)

    # Calculate properties
    properties_list = calculate(filename, scr=scr)

    return properties_list


def calculate_parameters(
    filename,
    parameters_list,
    scr=None,
    n_procs=1):
    """

    # TODO Create uuid folders in scr
    # TODO Change dir for each thread
    # TODO Stream results for each thread

    # TODO Collect properties, same order
    """

    if scr is None:
        scr = "./"

    Path(scr).mkdir(parents=True, exist_ok=True)

    proc_dirs = range(1, n_procs+1)
    proc_dirs = list(proc_dirs)
    proc_dirs = [os.path.join(scr, str(d)) for d in proc_dirs]

    # Prepare scrdirs with input
    for d in proc_dirs:
        Path(d).mkdir(parents=True, exist_ok=True)
        shutil.copy2(os.path.join(scr, filename), os.path.join(d, filename))

    kwargs = {
        "scr": scr,
        "filename": filename,
    }

    mapfunc = functools.partial(param_worker, **kwargs)

    p = mp.Pool(n_procs)
    results = p.map(mapfunc, parameters_list)

    return results


def numerical_jacobian(inputstr, params, dh=10**-5, n_procs=2):
    """
    get properties for

    """

    params_joblist = []
    atom_keys = params.keys()
    param_keys = {}
    param_grad = {}

    for atom in atom_keys:
        keys = params.get(atom).keys()
        param_keys[atom] = keys

        param_grad[atom] = {}
        for key in keys:
            param_grad[atom][key] = []

    for atom in atom_keys:
        for key in param_keys[atom]:

            dparams = copy.deepcopy(params)

            # forward
            dparams[atom][key] += dh
            params_joblist.append(copy.deepcopy(dparams))

            # backward
            dparams[atom][key] -= 2*dh
            params_joblist.append(copy.deepcopy(dparams))


    # Calculate all results
    results = calculate_multi_params(inputstr, params_joblist, n_procs=n_procs)
    n_results = len(results)

    i = 0
    for atom in atom_keys:
        for key in param_keys[atom]:
            param_grad[atom][key].append(results[i])
            param_grad[atom][key].append(results[i+1])
            i += 2

    return param_grad


def set_params(parameters, scr=None, **kwargs):
    """
    """

    txt = dump_params(parameters, **kwargs)

    filename = "fort.14"

    if scr is not None:
        filename = os.path.join(scr, filename)

    with open(filename, 'w') as f:
        f.write(txt)

    return


def load_params(filename, ignore_keys=IGNORE_PARAMS):

    with open(filename, 'r') as f:
        params = f.read()
        params = json.loads(params)

    if ignore_keys is not None:
        for atom in params:
            for key in ignore_keys:
                params[atom].pop(key, None)

    return params


def dump_params(parameters, ignore_keys=IGNORE_PARAMS):
    """
    """

    txt = ""

    for atomtype in parameters:
        for key in parameters[atomtype]:
            if key in ignore_keys:
                continue
            line = "{:8s} {:2s} {:15.11f}\n".format(key, atomtype, parameters[atomtype][key])
            txt += line

    return txt


def get_inputs(atoms_list, coords_list, charges, titles, **kwargs):
    """
    """

    inptxt = ""
    for atoms, coords, charge, title in zip(atoms_list, coords_list, charges, titles):
        txt = get_input(atoms, coords, charge, title=title, **kwargs)
        inptxt += txt

    return inptxt


def get_input(atoms, coords, charge,
    title="",
    method="PM3",
    header=DEFAULT_HEADER,
    read_params=True,
    **kwargs):
    """
    """

    # WARNING: INTERNAL COORDINATES ARE ASSUMED -
    # FOR THREE-ATOM SYSTEMS

    n_atoms = len(atoms)

    txt = header.format(method=method, charge=charge, title=title)

    if read_params:
        txt = txt.split("\n")
        txt[0] += " iparok=1"
        txt = "\n".join(txt)

    txt += "\n"

    if n_atoms <= 3:
        txt += internal_coordinates(atoms, coords, **kwargs)

    else:

        for atom, coord in zip(atoms, coords):
            line = ATOMLINE.format(atom, *coord)
            txt += line + "\n"

    txt += "\n"

    return txt


def internal_coordinates(atomtypes, xyz, optimize=False):

    natoms = len(atomtypes)

    opt_flag = 0
    if optimize: opt_flag = 1
    output = ""

    if (natoms == 3):

        ba = xyz[1] - xyz[0]
        bc = xyz[1] - xyz[2]
        cosine_angle = np.dot(ba, bc) / (np.linalg.norm(ba) * np.linalg.norm(bc))
        angle = np.arccos(cosine_angle) / np.pi * 180.0

        norm_ba = np.linalg.norm(ba)
        norm_bc = np.linalg.norm(bc)

        output += f"{atomtypes[0]}\n"
        output += f"{atomtypes[1]} {norm_ba} {opt_flag}\n"
        output += f"{atomtypes[2]} {norm_bc} {opt_flag} {angle} {opt_flag}\n"

    elif (natoms == 2):

        ba = xyz[1] - xyz[0]
        norm_ba = np.linalg.norm(ba)
        output += f"{atomtypes[0]}\n"
        output += f"{atomtypes[1]} {norm_ba} {opt_flag}\n"

    elif (natoms == 1):

        output += f"{atomtypes[0]}\n"

    return output


def get_default_params(method):
    """

    Get the default parameters of a method

    # TODO Extend the atom types

    """

    atoms = [x.strip().upper() for x in [
        'h ', \
        'c ', 'n ', 'o ', 'f ',\
        ]]

    n_atoms = len(atoms)
    method = method.upper()
    header = "{:} 0SCF MULLIK PRECISE charge={{:}} jprint=5\n\nTITLE {{:}}".format(method)

    coords = np.arange(n_atoms*3)
    coords = coords.reshape((n_atoms, 3))
    coords *= 5

    txt = get_input(atoms, coords, 0, "get params", header=header)
    filename = "_tmp_params.inp"

    with open(filename, 'w') as f:
        f.write(txt)

    molecules = run_mndo_file(filename)

    lines = next(molecules)

    idx = misc.get_index(lines, "PARAMETER VALUES USED IN THE CALCULATION")
    idx += 4

    parameters = {}

    while True:

        line = lines[idx]
        line = line.strip().split()

        if len(line) == 0:

            line = lines[idx+1]
            line = line.strip().split()

            if len(line) == 0:
                break
            else:
                idx += 1
                continue

        atom = line[0]
        param = line[1]
        value = line[2]
        unit = line[3]
        desc = " ".join(line[4:])

        atom = int(atom)
        atom = convert_atom(atom)
        value = float(value)

        if atom not in parameters.keys():
            parameters[atom] = {}

        parameters[atom][param] = value

        idx += 1

    return parameters



def dump_default_parameters():
    """
    Dump mndo.exe method parameters (json format) on disk

    Will be stored as parameters.method_name.json.

    """

    methods = ["MNDO", "AM1", "PM3", "OM2"]

    for method in methods:
        parameters = get_default_params(method)
        filename = "parameters.{:}.json".format(method.lower())
        with open(filename, 'w') as f:
            json.dump(parameters, f, indent=4)

    return


def write_input_file(
    atoms_list,
    coords_list,
    charges,
    titles,
    method,
    filename,
    **kwargs):

    txt = get_inputs(atoms_list, coords_list, charges, titles, method=method, **kwargs)

    with open(filename, "w") as f:
        f.write(txt)

    return


def main():

    # # Load data
    # atoms_list, coord_list, charges, titles = load_data()
    #
    # # TODO Select data
    #
    # # TODO Set input file
    # txt = get_inputs(atoms_list, coord_list, charges, titles)
    # f = open("runfile.inp", 'w')
    # f.write(txt)
    # f.close()

    # TODO set params
    # set_params(__PARAMETERS_OM2__)

    # TODO Run calculation
    # stdout = run_mndo_file("runfile.inp")
    #
    # for lines in stdout:
    #     properties = get_properties(lines)
    #     print(properties)

    # TODO Parse properties

    # Test multi input
    params = load_params("parameters/parameters.pm3.json")
    params.pop("F")

    # params_list = [params]*100

    with open("_tmp_optimizer") as f:
        inputstr = f.read()

    # results = calculate_multi_params(inputstr, params_list, n_procs=1)

    params_grad = numerical_jacobian(inputstr, params)
    print(json.dumps(params_grad, indent=1))

    return


if __name__ == '__main__':
    # Chemhelp
    main()


