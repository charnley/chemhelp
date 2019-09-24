"""

Module for chemical informatic tasks

"""

from io import StringIO
import sys
import gzip

import numpy as np

import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdkit.Chem.Draw as Draw
import rdkit.Chem.ChemicalForceFields as ChemicalForceFields
import rdkit.Chem.rdMolDescriptors as rdMolDescriptors

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


def molobj_to_atoms(molobj, atom_type="int"):

    atoms = molobj.GetAtoms()

    if atom_type == "str":
        atoms = [atom.GetSymbol() for atom in atoms]

    elif atom_type == "int":
        atoms = [atom.GetAtomicNum() for atom in atoms]
        atoms = np.array(atoms)

    return atoms


def molobj_to_xyz(molobj, atom_type="int"):
    """
    rdkit molobj to xyz
    """

    atoms = molobj_to_atoms(molobj, atom_type=atom_type)

    coordinates = molobj_to_coordinates(molobj)

    return atoms, coordinates


def molobj_add_hydrogens(molobj):

    molobj = Chem.AddHs(molobj)

    return molobj

def molobj_optimize(molobj):

    status = AllChem.EmbedMolecule(molobj)
    status = AllChem.UFFOptimizeMolecule(molobj)

    return status


def molobj_to_sdfstr(mol):
    """

    there must be a easier way to do this

    """

    # Chem rdkit::MolToMolBlock ???!
    # txt = Chem.MolToMolBlock(mol)

    sio = StringIO()
    w = Chem.SDWriter(sio)
    w.write(mol)
    w.flush()
    sdfstr = sio.getvalue()

    return sdfstr


def molobj_to_smiles(mol, remove_hs=False):
    """

    RDKit Mol Obj to SMILES wrapper

    """
    if remove_hs:
        mol = Chem.RemoveHs(mol)

    smiles = Chem.MolToSmiles(mol)

    return smiles


def molobj_to_svgstr(molobj,
                     highlights=None,
                     pretty=False,
                     removeHs=False):
    """

    Returns SVG in string format

    """

    if removeHs:
        molobj = Chem.RemoveHs(molobj)

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


def sdfstr_to_molobj(sdfstr, remove_hs=False):
    """
    SDF to mol obj
    """

    # sio = sys.stderr = StringIO()
    mol = Chem.MolFromMolBlock(sdfstr, removeHs=remove_hs)

    # if mol is None:
    #     return None, sio.getvalue()

    return mol, ""


def sdfstr_to_smiles(sdfstr, remove_hs=False):
    """
    SDF to SMILES converter
    """

    # sio = sys.stderr = StringIO()
    mol = Chem.MolFromMolBlock(sdfstr, removeHs=remove_hs)

    if mol is None:
        return None, "" #sio.getvalue()

    smiles = Chem.MolToSmiles(mol)
    status = ""

    return smiles, status


def smiles_to_sdfstr(smilesstr, add_hydrogens=True):
    """
    SMILES to SDF converter
    """

    # sio = sys.stderr = StringIO()
    mol = Chem.MolFromSmiles(smilesstr)

    if mol is None:
        return None, "" #sio.getvalue()

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    sdfstr = molobj_to_sdfstr(mol)
    status = ""

    return sdfstr, status


def smiles_to_molobj(smilesstr, add_hydrogens=True):

    mol = Chem.MolFromSmiles(smilesstr)

    if add_hydrogens:
        mol = Chem.AddHs(mol)

    return mol, ""


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


def genereate_conformers(smilesstr, max_conf=20, min_conf=10):

    molobj, status = smiles_to_molobj(smilesstr, add_hydrogens=True)

    if molobj is None:
        return None

    status = AllChem.EmbedMolecule(molobj)
    status = AllChem.UFFOptimizeMolecule(molobj)

    rot_bond = rdMolDescriptors.CalcNumRotatableBonds(molobj)

    confs = min(1 + 3*rot_bond, max_conf)
    confs = max(confs, min_conf)

    AllChem.EmbedMultipleConfs(molobj, numConfs=confs,
                useExpTorsionAnglePrefs=True,
                useBasicKnowledge=True)

    return molobj


def conformationalsearch(smiles):

    molobj = genereate_conformers(smiles)

    if molobj is None:
        return None

    conformers = molobj.GetConformers()

    # status will be 0 for converged molecules
    res = AllChem.MMFFOptimizeMoleculeConfs(molobj)
    res = np.array(res)

    status = res[:,0]
    energies = res[:,1]
    idx_converged, = np.where(status == 0)
    energies = energies[idx_converged]

    if energies.shape[0] == 0:
        return None

    idx_lowest = np.argsort(energies)[0]
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


if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('--sdf', action='store', help='', metavar='FILE')
    parser.add_argument('--smi', action='store', help='', metavar='FILE')

    args = parser.parse_args()

    if args.smi is not None:
        molobj = conformationalsearch(args.smi)
        sdfstr = molobj_to_sdfstr(molobj)
        print(sdfstr)


    pass

