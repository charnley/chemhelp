"""

Module for chemical informatic tasks

TODO AllChem.Compute2DCoords(m) for png

"""

import os
from io import StringIO
import sys
import gzip

import numpy as np

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors
import rdkit.Chem.rdmolops as rdmolops

# Chem.WrapLogs()

# spin-multiplicities 2,3,4,3,2 for the atoms H, C, N, O, F, respectively.
MULTIPLICITY = {}
MULTIPLICITY["H"] = 2
MULTIPLICITY["C"] = 3
MULTIPLICITY["N"] = 4
MULTIPLICITY["O"] = 3
MULTIPLICITY["F"] = 2
MULTIPLICITY["Cl"] = 2

ATOM_LIST = [x.strip() for x in [
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
        idx = ATOM_LIST.index(atom) + 1
        return idx

    else:
        atom = ATOM_LIST[atom -1].capitalize()
        return atom


def atom_str(iatm):

    atom = ATOM_LIST[iatm-1]
    atom = atom.capitalize()

    return atom


def atom_int(atmstr):
    atom = atmstr.lower()
    atom = ATOM_LIST.index(atom) + 1
    return atom


def clean_sdf_header(sdfstr):

    sdfstr = str(sdfstr)
    for _ in range(2):
        i = sdfstr.index('\n')
        sdfstr = sdfstr[i+1:]
    sdfstr = "\n\n" + sdfstr

    return sdfstr


def read_sdffile(filename, remove_hs=False, sanitize=True):
    """
    """

    ext = filename.split(".")[-1]

    if ext == "sdf":

        suppl = Chem.SDMolSupplier(filename,
            removeHs=remove_hs,
            sanitize=sanitize)

    elif ext == "gz":

        fobj = gzip.open(filename)
        suppl = Chem.ForwardSDMolSupplier(fobj,
            removeHs=remove_hs,
            sanitize=sanitize)

    else:
        print("could not read file")
        quit()

    return suppl


def read_xyzfile(filename):


    return


def get_torsions(mol):
    """ return idx of all torsion pairs
    All heavy atoms, and one end can be a hydrogen
    """

    any_atom = "[*]"
    not_hydrogen = "[!H]"

    smarts = [
        any_atom,
        any_atom,
        any_atom,
        any_atom]

    smarts = "~".join(smarts)

    idxs = mol.GetSubstructMatches(Chem.MolFromSmarts(smarts))
    idxs = [list(x) for x in idxs]
    idxs = np.array(idxs)

    rtnidxs = []

    for idx in idxs:

        atoms = get_torsion_atoms(mol, idx)
        atoms = np.array(atoms)
        idxh, = np.where(atoms == "H")

        if idxh.shape[0] > 1: continue
        elif idxh.shape[0] > 0:
            if idxh[0] == 1: continue
            if idxh[0] == 2: continue

        rtnidxs.append(idx)

    return np.array(rtnidxs, dtype=int)


def get_torsion_atoms(mol, torsion):
    """
    return all atoms for specific torsion indexes

    # TODO Not really need if you think it through

    """

    atoms = mol.GetAtoms()
    atoms = [atom.GetSymbol() for atom in atoms]
    atoms = np.array(atoms)
    atoms = atoms[torsion]

    return atoms


def canonical(smiles):
    """
    Translate smiles into a canonical form
    """

    molobj = Chem.MolFromSmiles(smiles)
    smiles = Chem.MolToSmiles(molobj, canonical=True)

    return smiles


def molobj_to_coordinates(molobj):

    conformer = molobj.GetConformer()
    coordinates = conformer.GetPositions()
    coordinates = np.array(coordinates)

    return coordinates


def molobj_to_atoms(molobj, atom_type=int):

    atoms = molobj.GetAtoms()

    if atom_type == str or atom_type == "str":
        atoms = [atom.GetSymbol() for atom in atoms]

    elif atom_type == int or atom_type == "int":
        atoms = [atom.GetAtomicNum() for atom in atoms]
        atoms = np.array(atoms)

    return atoms


def molobj_to_xyz(molobj, atom_type=int):
    """
    rdkit molobj to xyz
    """

    atoms = molobj_to_atoms(molobj, atom_type=atom_type)

    coordinates = molobj_to_coordinates(molobj)

    return atoms, coordinates


def molobj_to_axyzc(molobj, atom_type="int"):
    """
    rdkit molobj to atoms, xyz, charge
    """

    atoms = molobj_to_atoms(molobj, atom_type=atom_type)

    coordinates = molobj_to_coordinates(molobj)

    charge = rdmolops.GetFormalCharge(molobj)

    return atoms, coordinates, charge


def molobj_add_hydrogens(molobj):

    molobj = Chem.AddHs(molobj)

    return molobj

def molobj_optimize(molobj, max_steps=1000):

    status_embed = AllChem.EmbedMolecule(molobj)

    if status_embed != 0:
        return status_embed

    try:
        status_2 = AllChem.UFFOptimizeMolecule(molobj, maxIters=max_steps)
    except RuntimeError:
        status_2 = 5

    return status_2


def molobj_to_sdfstr(mol, use_v3000=False):
    """

    there must be a easier way to do this

    """

    if mol is None: return None

    n_confs = mol.GetNumConformers()

    txts = []

    for i in range(n_confs):
        txt = Chem.MolToMolBlock(mol, confId=i, forceV3000=use_v3000)
        txts += [txt]

    txts = "$$$$\n".join(txts)

    return txts


def molobj_to_smiles(mol, remove_hs=False):
    """

    RDKit Mol Obj to SMILES wrapper

    """
    if remove_hs:
        mol = Chem.RemoveHs(mol)

    smiles = Chem.MolToSmiles(mol)

    return smiles


def molobj_to_svgstr(molobj,
    force2d=False,
    highlights=None,
    pretty=False,
    removeHs=False):
    """

    Returns SVG in string format

    """

    if removeHs:
        molobj = Chem.RemoveHs(molobj)

    if force2d:
        molobj = copy.deepcopy(molobj)
        # TODO complete this func

    svg = Draw.MolsToGridImage(
        [molobj],
        molsPerRow=1,
        subImgSize=(400,400),
        useSVG=True,
        highlightAtomLists=[highlights])

    svg = svg.replace("xmlns:svg", "xmlns")

    if pretty:

        svg = svg.split("\n")

        for i, line in enumerate(svg):

            # Atom letters
            if "text" in line:

                replacetext = "font-size"
                borderline = "fill:none;fill-opacity:1;stroke:#FFFFFF;stroke-width:10px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;"

                # Add border to text
                border_text = line
                border_text = border_text.replace('stroke:none;', '')
                border_text = border_text.replace(replacetext, borderline+replacetext )

                svg[i] = border_text + "\n" + line

                continue


            if "path" in line:

                # thicker lines
                line = line.replace('stroke-width:2px', 'stroke-width:3px')
                svg[i] = line

        svg = "\n".join(svg)

    return svg


def sdfstr_to_molobj(sdfstr, remove_hs=False, return_status=False):
    """
    SDF to mol obj
    """

    if return_status:
        Chem.WrapLogs()
        sio = sys.stderr = StringIO()

    mol = Chem.MolFromMolBlock(sdfstr, removeHs=remove_hs)

    if return_status:
        msg = sio.getvalue()
        return mol, msg

    return mol


def sdfstr_to_smiles(sdfstr, remove_hs=False):
    """
    SDF to SMILES converter
    """

    mol = Chem.MolFromMolBlock(sdfstr, removeHs=remove_hs)

    smiles = Chem.MolToSmiles(mol)

    return smiles


def smiles_to_sdfstr(smilesstr,
    add_hydrogens=True,
    return_status=False):
    """
    SMILES to SDF converter
    """

    if return_status:
        Chem.WrapLogs()
        sio = sys.stderr = StringIO()

    mol = Chem.MolFromSmiles(smilesstr)

    if return_status:
        msg = sio.getvalue()

    if mol is not None and add_hydrogens:
        mol = Chem.AddHs(mol)

    sdfstr = molobj_to_sdfstr(mol)

    if return_status:
        sdfstr, msg

    return sdfstr


def smiles_to_molobj(smilesstr,
    add_hydrogens=True,
    return_status=False):

    if return_status:
        Chem.WrapLogs()
        sio = sys.stderr = StringIO()

    mol = Chem.MolFromSmiles(smilesstr)

    if return_status:
        msg = sio.getvalue()

    if mol is not None and add_hydrogens:
        mol = Chem.AddHs(mol)

    if return_status:
        return mol, msg

    return mol


def add_conformer(molobj, coordinates):

    conf = Chem.Conformer(len(coordinates))

    for i, coordinate in enumerate(coordinates):
        conf.SetAtomPosition(i, coordinate)

    molobj.AddConformer(conf, assignId=True)

    return


def molobj_copy(molobj):

    molobj_prime = Chem.Mol(molobj)

    return molobj_prime


def molobj_get_coordinates(molobj):
    """
    """

    conformer = molobj.GetConformer()
    coordinates = conformer.GetPositions()
    coordinates = np.asarray(coordinates)

    return coordinates


def conformer_set_coordinates(conformer, coordinates):

    for i, pos in enumerate(coordinates):
        conformer.SetAtomPosition(i, pos)

    return


def molobj_set_coordinates(molobj, coordinates):

    conformer = molobj.GetConformer()
    conformer_set_coordinates(conformer, coordinates)

    return


def generate_conformers(smilesstr, max_conf=20, min_conf=10, max_steps=1000):

    molobj = smiles_to_molobj(smilesstr, add_hydrogens=True)

    if molobj is None:
        return None

    status_embed = AllChem.EmbedMolecule(molobj)

    if status_embed != 0:
        return None

    status_optim = AllChem.UFFOptimizeMolecule(molobj, maxIters=max_steps)

    # Keep unconverged uff
    # if status_optim != 0:
    #     return None

    # Check bond lengths
    dist = Chem.rdmolops.Get3DDistanceMatrix(molobj)
    np.fill_diagonal(dist, 10.0)
    min_dist = np.min(dist)

    # For some atom_types in UFF, it will fail
    if min_dist < 0.001:
        print("fail", smilesstr)
        return None

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(molobj)

    confs = min(1 + 3*rot_bond, max_conf)
    confs = max(confs, min_conf)

    status = AllChem.EmbedMultipleConfs(molobj,
        numConfs=confs,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True)

    return molobj


def check_conformer_dist(molobj, cutoff=0.001):
    """
    For some atom_types in UFF, rdkit will fail optimization and stick multiple atoms ontop of eachother

    Especially in CS(F3)

    """

    n_confs = molobj.GetNumConformers()

    status = []

    for i in range(n_confs):
        dist = Chem.rdmolops.Get3DDistanceMatrix(molobj, confId=i)
        np.fill_diagonal(dist, 10.0)
        min_dist = np.min(dist)

        this = 0
        if min_dist < cutoff:
            this += 1

        status.append(this)

    status = np.array(status)

    return status


def conformationalsearch(smiles):

    molobj = genereate_conformers(smiles)

    if molobj is None:
        return None

    conformers = molobj.GetConformers()

    # status will be 0 for converged molecules
    res = AllChem.UFFOptimizeMoleculeConfs(molobj)
    res = np.array(res)

    try:
        status = res[:,0]
    except IndexError:
        return None

    # check distances
    status_dist = check_conformer_dist(molobj)
    status += status_dist

    idx_converged, = np.where(status == 0)

    energies = res[:,1]
    energies = energies[idx_converged]

    if energies.shape[0] == 0:
        return None

    # lowest converged energy
    idx_lowest = np.argsort(energies)[0]

    # convert to global idx
    idx_lowest = idx_converged[idx_lowest]

    coord = conformers[idx_lowest].GetPositions()

    molobj.RemoveAllConformers()
    add_conformer(molobj, coord)

    return molobj


def save_molobj(molobj, coordinates=None):
    """
    save sdf from specific coordinates
    TODO Move to molobj_to_sdfstr, optional
    """

    if coordinates is not None:
        conformer = molobj.GetConformer()
        conformer_set_coordinates(conformer, coordinates)

    sdf = molobj_to_sdfstr(molobj)

    return sdf


def molobj_to_mol2(molobj, charges=None):
    """
    function from
    https://www.mdanalysis.org/docs/_modules/MDAnalysis/coordinates/MOL2.html

    """

    # Bonds
    bond_lines = ["@<TRIPOS>BOND"]
    bond_fmt = "{0:>5} {1:>5} {2:>5} {3:>2}"
    bonds = list(molobj.GetBonds())
    n_bonds = len(bonds)
    for i, bond in enumerate(bonds):
        a = bond.GetBeginAtomIdx()
        b = bond.GetEndAtomIdx()
        t = bond.GetBondType()
        tf = bond.GetBondTypeAsDouble()

        if tf.is_integer():
            t = int(t)
        else:
            t = "ar"

        bond = bond_fmt.format(i+1, a+1, b+1, t)
        bond_lines.append(bond)

    bond_lines.append("\n")
    bond_lines = "\n".join(bond_lines)

    # Atoms
    atom_lines = ["@<TRIPOS>ATOM"]
    atom_fmt = "{0:>4} {1:>4} {2:>13.4f} {3:>9.4f} {4:>9.4f} {5:>4} {6} {7} {8:>7.4f}"
    atoms = list(molobj.GetAtoms())
    atoms_int = [atom.GetAtomicNum() for atom in atoms]
    atoms_str = [atom.GetSymbol() for atom in atoms]
    atoms_int = np.array(atoms_int)
    n_atoms = len(atoms)
    conformer = molobj.GetConformer()
    coordinates = conformer.GetPositions()
    coordinates = np.array(coordinates)
    unique_atoms = np.unique(atoms_int)

    if charges is None:
        charges = np.zeros(n_atoms)

    atm_i = 1


    for j in range(n_atoms):

        name = atoms_str[j]
        pos0 = coordinates[j,0]
        pos1 = coordinates[j,1]
        pos2 = coordinates[j,2]
        typ = atoms_str[j]
        resid = 0
        resname = "MOL"
        charge = charges[j]

        atmstr = atom_fmt.format(j+1, name, pos0, pos1, pos2, typ, resid,  resname, charge)
        atom_lines.append(atmstr)

        atm_i += 1
        continue

    atom_lines.append("")
    atom_lines = "\n".join(atom_lines)

    # Complete
    checksumstr = f"{n_atoms} {n_bonds} 0 0 0"
    head_lines = ["@<TRIPOS>MOLECULE", "TITLE"]
    head_lines += [checksumstr, "SMALL", "MULLIKEN_CHARGES", "NAME"]
    head_lines.append("")
    head_lines = "\n".join(head_lines)

    rtnstr = head_lines + atom_lines + bond_lines

    return rtnstr


def remove_salt(smiles):

    return max(smiles.split("."), key=len)


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('--smi', action='store', help='', metavar='FILE')

    args = parser.parse_args()

    if args.smi is not None:

        if os.path.exists(args.smi):

            with open(args.smi, 'r') as f:
                lines = f.readlines()
                lines = [line.strip() for line in lines]

            rtnstr = ""
            for line in lines:
                molobj = conformationalsearch(line)
                sdfstr = molobj_to_sdfstr(molobj)
                rtnstr += sdfstr
            print(rtnstr)

        else:

            molobj = conformationalsearch(args.smi)
            sdfstr = molobj_to_sdfstr(molobj)
            print(sdfstr)


    pass

