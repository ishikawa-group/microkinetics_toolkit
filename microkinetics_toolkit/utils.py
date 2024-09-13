import numpy as np
from ase import Atom, Atoms


def read_reactionfile(file):
    """
    Read reaction file and return reactants, reactions, and products.

    Args:
        file (str): reaction file name
    Returns:
        reac (list): reactants
        rxn (list): reaction type
        prod (list): products
    """
    import re

    f = open(file, "r")

    # drop comment and branck lines
    lines = f.readlines()
    newlines = []
    for line in lines:
        if not (re.match(r"^#", line)) and not (re.match(r"^s*$", line)):
            newlines.append(line)
    lines = newlines
    numlines = len(lines)

    reac = list(range(numlines))
    rxn = list(range(numlines))
    prod = list(range(numlines))

    for i, line in enumerate(lines):
        text = line.replace("\n", "").replace(">", "").split("--")
        reac_tmp = text[0]
        rxn_tmp = text[1]
        prod_tmp = text[2]

        reac[i] = re.split(" \+ ", reac_tmp)  # for cations
        prod[i] = re.split(" \+ ", prod_tmp)  # for cations

        reac[i] = remove_space(reac[i])
        prod[i] = remove_space(prod[i])

        rxn[i] = reac[i][0] + "_" + rxn_tmp

    return reac, rxn, prod


def return_lines_of_reactionfile(file):
    """
    Return lines of reaction file.
    """
    import re

    # drop comment and branck lines
    f = open(file, "r")
    lines = f.readlines()
    newlines = []
    for line in lines:
        if not (re.match(r"^#", line)) and not (re.match(r"^s*$", line)):
            newlines.append(line)
    lines = newlines
    return lines


def remove_space(obj):
    """
    Remove space in the object.
    """
    newobj = [0] * len(obj)
    if isinstance(obj, str):
        #
        # string
        #
        newobj = obj.replace(" ", "")
    elif isinstance(obj, list):
        #
        # list
        #
        for i, obj2 in enumerate(obj):
            if isinstance(obj2, list):
                #
                # nested list
                #
                for ii, jj in enumerate(obj2):
                    jj = jj.strip()
                newobj[i] = jj
            elif isinstance(obj2, str):
                #
                # simple list
                #
                obj2 = obj2.replace(" ", "")
                newobj[i] = obj2
            elif isinstance(obj2, int):
                #
                # integer
                #
                newobj[i] = obj2
            else:
                newobj[i] = obj2
    else:  # error
        print("remove_space: input str or list")

    return newobj


def get_reac_and_prod(reactionfile):
    """
    Form reactant and product information.
    """
    (reac, rxn, prod) = read_reactionfile(reactionfile)

    rxn_num = len(rxn)

    r_ads = list(range(rxn_num))
    r_site = [[] for _ in range(rxn_num)]
    r_coef = [[] for _ in range(rxn_num)]

    p_ads = list(range(rxn_num))
    p_site = list(range(rxn_num))
    p_coef = list(range(rxn_num))

    for irxn, jrnx in enumerate(rxn):
        ireac = reac[irxn]
        iprod = prod[irxn]
        ireac_num = len(ireac)
        iprod_num = len(iprod)
        #
        # reactant
        #
        r_ads[irxn] = list(range(ireac_num))
        r_site[irxn] = list(range(ireac_num))
        r_coef[irxn] = list(range(ireac_num))

        for imol, mol in enumerate(ireac):
            r_site[irxn][imol] = []
            r_ads[irxn][imol] = []
            #
            # coefficient
            #
            if "*" in mol:
                # r_coef[irxn][imol] = int(mol.split("*")[0])
                r_coef[irxn][imol] = float(mol.split("*")[0])
                rest = mol.split("*")[1]
            else:
                r_coef[irxn][imol] = 1
                rest = mol

            # site
            if ',' in rest:
                sites = rest.split(',')
                for isite, site in enumerate(sites):
                    r_site[irxn][imol].append(site.split('_')[1])
                    r_ads[irxn][imol].append(site.split('_')[0])
            elif '_' in rest:
                r_site[irxn][imol].append(rest.split('_')[1])
                r_ads[irxn][imol].append(rest.split('_')[0])
            else:
                r_site[irxn][imol].append('gas')
                r_ads[irxn][imol].append(rest)
        #
        # product
        #
        p_ads[irxn] = list(range(iprod_num))
        p_site[irxn] = list(range(iprod_num))
        p_coef[irxn] = list(range(iprod_num))

        for imol, mol in enumerate(iprod):
            p_site[irxn][imol] = []
            p_ads[irxn][imol] = []
            #
            # coefficient
            #
            if "*" in mol:
                # p_coef[irxn][imol] = int(mol.split("*")[0])
                p_coef[irxn][imol] = float(mol.split("*")[0])
                rest = mol.split("*")[1]
            else:
                p_coef[irxn][imol] = 1
                rest = mol

            # site
            if ',' in rest:
                sites = rest.split(',')
                for isite, site in enumerate(sites):
                    p_site[irxn][imol].append(site.split('_')[1])
                    p_ads[irxn][imol].append(site.split('_')[0])
            elif '_' in rest:
                p_site[irxn][imol].append(rest.split('_')[1])
                p_ads[irxn][imol].append(rest.split('_')[0])
            else:
                p_site[irxn][imol].append('gas')
                p_ads[irxn][imol].append(rest)

    return r_ads, r_site, r_coef, p_ads, p_site, p_coef


def get_number_of_reaction(reactionfile):
    """
    Return number of elementary reactions
    """
    (reac, rxn, prod) = read_reactionfile(reactionfile)
    rxn_num = len(rxn)
    return rxn_num


def get_preexponential(reactionfile):
    """
    Return pre-exponential factor.
    Not needed for MATLAB use
    """
    from ase.collections import methane

    #
    # calculate pre-exponential factor
    #
    (r_ads, r_site, r_coef, p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

    rxn_num = get_number_of_reaction(reactionfile)

    Afor = np.array(rxn_num, dtype="f")
    Arev = np.array(rxn_num, dtype="f")

    mass_reac = np.array(rxn_num * [range(len(r_ads[0]))], dtype="f")
    mass_prod = np.array(rxn_num * [range(len(p_ads[0]))], dtype="f")

    for irxn in range(rxn_num):
        #
        # reactants
        #
        for imol, mol in enumerate(r_ads[irxn]):
            tmp = methane[mol]
            mass = sum(tmp.get_masses())
            mass_reac[irxn, imol] = mass
        #
        # products
        #
        for imol, mol in enumerate(p_ads[irxn]):
            tmp = methane[mol]
            mass = sum(tmp.get_masses())
            mass_prod[irxn, imol] = mass

        Afor[irxn] = 1.0
        Arev[irxn] = 1.0

    return Afor, Arev


def get_rate_coefficient(reactionfile, Afor, Arev, Efor, Erev, T):
    """
    Return rate coefficient.
    Not needed for MATLAB use
    """
    # calculate rate constant
    # (r_ads, r_site, r_coef,  p_ads, p_site, p_coef) = get_reac_and_prod(reactionfile)

    rxn_num = get_number_of_reaction(reactionfile)

    kfor = np.array(rxn_num, dtype="f")
    krev = np.array(rxn_num, dtype="f")

    RT = 8.314 * T / 1000.0  # in kJ/mol

    for irxn in range(rxn_num):
        kfor[irxn] = Afor[irxn] * np.exp(-Efor[irxn] / RT)
        krev[irxn] = Arev[irxn] * np.exp(-Erev[irxn] / RT)

    return kfor, krev


def read_speciesfile(speciesfile):
    """
    read species
    """
    f = open(speciesfile)

    species = f.read()
    species = species.replace('[', '')
    species = species.replace(']', '')
    species = species.replace(' ', '')
    species = species.replace('\n', '')
    species = species.replace('\'', '')
    species = species.split(",")

    return species


def remove_parentheses(file):
    """
    remove parentheses -- maybe for MATLAB use
    """
    import os
    tmpfile = "ttt.txt"
    os.system('cat %s | sed "s/\[//g" > %s' % (file, tmpfile))
    os.system('cat %s | sed "s/\]//g" > %s' % (tmpfile, file))
    os.system('rm %s' % tmpfile)


def get_species_num(*species):
    """
    Return what is the number of species in speciesfile.
    If argument is not present, returns the number of species.
    """
    from reaction_tools import read_speciesfile

    speciesfile = "species.txt"
    lst = read_speciesfile(speciesfile)

    if len(species) == 0:
        # null argument: number of species
        return len(lst)
    else:
        # return species number
        spec = species[0]
        return lst.index(spec)


def get_adsorption_sites(infile):
    """
    Read adsorption sites.
    """
    from reaction_tools import remove_space

    f = open(infile, "r")

    lines = f.readlines()

    mol = list(range(len(lines)))
    site = list(range(len(lines)))

    for i, line in enumerate(lines):
        aaa, bbb = line.replace("\n", "").split(":")
        mol[i] = remove_space(aaa)
        bbb = remove_space(bbb)
        site[i] = bbb.split(",")

    return mol, site


def find_closest_atom(surf, offset=(0, 0)):
    """
    Find the closest atom to the adsorbate.
    """
    from ase.build import add_adsorbate

    dummy = Atom('H', (0, 0, 0))
    ads_height = 0.1
    add_adsorbate(surf, dummy, ads_height, position=(0, 0), offset=offset)
    natoms = len(surf.get_atomic_numbers())
    last = natoms - 1
    # ads_pos = surf.get_positions(last)

    dist = surf.get_distances(last, [range(natoms)], vector=False)
    dist = np.array(dist)
    dist = np.delete(dist, last)  # delete adsorbate itself
    clothest = np.argmin(dist)

    return clothest


def sort_atoms_by_z(atoms):
    """
    Sort atoms by z-coordinate.
    """
    import collections

    #
    # keep information for original Atoms
    #
    tags = atoms.get_tags()
    pbc = atoms.get_pbc()
    cell = atoms.get_cell()

    dtype = [("idx", int), ("z", float)]
    #
    # get set of chemical symbols
    #
    symbols = atoms.get_chemical_symbols()
    elements = sorted(set(symbols), key=symbols.index)
    num_elem = []
    for i in elements:
        num_elem.append(symbols.count(i))

    #
    # loop over each groups
    #
    iatm = 0
    newatoms = Atoms()
    zcount = []
    for inum in num_elem:
        zlist = np.array([], dtype=dtype)
        for idx in range(inum):
            tmp = np.array([(iatm, atoms[iatm].z)], dtype=dtype)
            zlist = np.append(zlist, tmp)
            iatm = iatm + 1

        zlist = np.sort(zlist, order="z")

        for i in zlist:
            idx = i[0]
            newatoms.append(atoms[idx])

        tmp = np.array([], dtype=float)
        for val in zlist:
            tmp = np.append(tmp, round(val[1], 2))
        tmp = collections.Counter(tmp)
        zcount.append(list(tmp.values()))
    #
    # restore tag, pbc, cell
    #
    newatoms.set_tags(tags)
    newatoms.set_pbc(pbc)
    newatoms.set_cell(cell)

    return newatoms, zcount


def get_number_of_valence_electrons(atoms):
    """
    Returns number of valence electrons for VASP calculation.
    Calls VASP calculation once.
    Returned electron numbers should be ++1 or --1 for cations, anions, etc.
    """
    from ase.calculators.vasp import Vasp
    atoms.calc = Vasp(prec="normal", xc="PBE", encut=300.0, nsw=0, nelm=1)
    nelec = atoms.calc.get_number_of_electrons()
    return nelec


def read_charge(mol):
    """
    Read charge from molecule.
    """
    charge = 0.0  # initial
    if "^" in mol:
        neutral = False
        mol, charge = mol.split('^')
        charge = float(charge.replace('{', '').replace('}', '').replace('+', ''))
    else:
        neutral = True

    return mol, neutral, charge


def remove_side_and_flip(mol):
    """
    Remove SIDE and FLIP in molecule
    """
    if '-SIDEx' in mol:
        mol = mol.replace('-SIDEx', '')
    elif '-SIDEy' in mol:
        mol = mol.replace('-SIDEy', '')
    elif '-SIDE' in mol:
        mol = mol.replace('-SIDE', '')
    elif '-FLIP' in mol:
        mol = mol.replace('-FLIP', '')
    elif '-TILT' in mol:
        mol = mol.replace('-TILT', '')
    elif '-HIGH' in mol:
        mol = mol.replace('-HIGH', '')

    return mol


def neb_copy_contcar_to_poscar(nimages):
    """
    Copy 0X/CONTCAR to 0X/POSCAR after NEB run.
    """
    import os
    for images in range(nimages):
        os.system('cp %02d/CONTCAR %02d/POSCAR' % (images + 1, images + 1))


def make_it_closer_by_exchange(atom1, atom2, thre=0.05):
    """
    Exchange atoms to make it closer.

    Args:
        atom1 (Atoms): Atoms object
        atom2 (Atoms): Atoms object
        thre: when distance is larger than this value, do switch
    """
    from ase.geometry import distance

    natoms = len(atom1)
    const_list = atom1.constraints[0].get_indices()
    # Constrained atoms. Do not exchange these.

    for i in range(natoms):
        for j in range(i + 1, natoms):
            if atom1[i].symbol == atom1[j].symbol:
                # measure distance between "ATOMS" (in ASE object)
                dis_bef = distance(atom1, atom2)
                atom1B = atom1.copy()
                atom1B.positions[[i, j]] = atom1B.positions[[j, i]]
                dis_aft = distance(atom1B, atom2)

                if (dis_aft + thre < dis_bef) and not (i in const_list) and not (j in const_list):
                    #
                    # switch
                    #
                    print("exchanged {0:3d} and {1:3d}: (dis_bef, dis_aft) = ({2:6.2f},{3:6.2f})".
                          format(i, j, dis_bef, dis_aft))
                    tmp = atom1B.numbers[i]
                    atom1B.numbers[i] = atom1B.numbers[j]
                    atom1B.numbers[j] = tmp

                    atom1 = atom1B.copy()

    return atom1


def get_adsorbate_type(adsorbate, site):
    """
    Returns adsorbate type.

    Args:
        adsorbate (str): Molecule name.
        site (str): Site name.
    Returns:
        adsorbate_type (str): Adsorbate type
            "gaseous": gaseous molecule
            "surface": bare surface
            "adsorbed": adsorbed molecule
    """
    if site == "gas":
        if "surf" in adsorbate:
            ads_type = "surface"
        else:
            ads_type = "gaseous"
    else:
        ads_type = "adsorbed"

    return ads_type


def make_surface_from_cif(cif_file, indices=[1, 0, 0], repeat=[1, 1, 1], vacuum=10.0):
    """
    Make a surface from a CIF file.
    """
    from ase.build import surface
    from ase.io import read

    bulk = read(cif_file)
    bulk = bulk*repeat
    surf = surface(bulk, indices=indices, layers=2, vacuum=vacuum)

    surf.translate([0, 0, -vacuum+0.1])
    surf.pbc = True

    return surf


def replace_element(atoms, from_element, to_element, percent=100):
    import random

    from ase.build import sort

    elements = atoms.get_chemical_symbols()
    num_from_elements = elements.count(from_element)
    num_replace = int((percent/100) * num_from_elements)

    indices = [i for i, j in enumerate(elements) if j == from_element]
    random_item = random.sample(indices, num_replace)
    for i in random_item:
        atoms[i].symbol = to_element

    atoms = sort(atoms)
    return atoms


def run_packmol(xyz_file, a, num, outfile):
    import os

    packmol = "/Users/ishi/packmol/packmol"
    filetype = "xyz"

    cell1 = [0.0, 0.0, 0.0, a, a, a]
    cell2 = " ".join(map(str, cell1))

    f = open("pack_tmp.inp", "w")
    text = ["tolerance 2.0"             + "\n",
            "output "     + outfile     + "\n",
            "filetype "   + filetype    + "\n",
            "structure "  + xyz_file    + "\n",
            "  number "   + str(num)    + "\n",
            "  inside box " + cell2     + "\n",
            "end structure"]
    f.writelines(text)
    f.close()

    run_string = packmol + " < pack_tmp.inp"

    os.system(run_string)

    # os.system("rm pack_tmp.inp")


def json_to_csv(jsonfile, csvfile):
    import json

    import pandas as pd
    from pandas.io.json import json_normalize
    f = open(jsonfile, "r")
    d = json.load(f)

    dd = []
    nrec = len(d)
    for i in range(1, nrec):
        if str(i) in d:
            tmp = d[str(i)]
            dd.append(json_normalize(tmp))

    ddd = pd.concat(dd)

    newcol = []
    for key in ddd.columns:
        key = key.replace("calculator_parameters.", "")
        key = key.replace("key_value_pairs.", "")
        key = key.replace("data.", "")
        newcol.append(key)

    ddd.columns = newcol

    # sort data by "num"
    if "num" in ddd.columns:
        ddd2 = ddd.set_index("num")
        ddd  = ddd2.sort_index()

    ddd.to_csv(csvfile)


def load_ase_json(jsonfile):
    import json

    import pandas as pd
    f = open(jsonfile, "r")
    d = json.load(f)

    dd = []
    nrec = len(d)
    for i in range(1, nrec):
        if str(i) in d:
            tmp = d[str(i)]
            dd.append(pd.json_normalize(tmp))

    ddd = pd.concat(dd)

    newcol = []
    for key in ddd.columns:
        key = key.replace("calculator_parameters.", "")
        key = key.replace("key_value_pairs.", "")
        key = key.replace("data.", "")
        newcol.append(key)

    ddd.columns = newcol

    # sort data by "num"
    if "num" in ddd.columns:
        ddd2 = ddd.set_index("num")
        ddd  = ddd2.sort_index()

    return ddd


def delete_num_from_json(num, jsonfile):
    from ase.db import connect

    db = connect(jsonfile)
    id_ = db.get(num=num).id
    db.delete([id_])


def sort_atoms_by(atoms, xyz="x"):
    # keep information for original Atoms
    tags = atoms.get_tags()
    pbc  = atoms.get_pbc()
    cell = atoms.get_cell()
    dtype = [("idx", int), (xyz, float)]

    newatoms = Atoms()
    symbols = list(set(atoms.get_chemical_symbols()))
    for symbol in symbols:
        subatoms = Atoms(list(filter(lambda x: x.symbol == symbol, atoms)))
        atomlist = np.array([], dtype=dtype)
        for idx, atom in enumerate(subatoms):
            if xyz == "x":
                tmp = np.array([(idx, atom.x)], dtype=dtype)
            elif xyz == "y":
                tmp = np.array([(idx, atom.y)], dtype=dtype)
            else:
                tmp = np.array([(idx, atom.z)], dtype=dtype)

            atomlist = np.append(atomlist, tmp)

        atomlist = np.sort(atomlist, order=xyz)

        for i in atomlist:
            idx = i[0]
            newatoms.append(subatoms[idx])

    # restore
    newatoms.set_tags(tags)
    newatoms.set_pbc(pbc)
    newatoms.set_cell(cell)

    return newatoms


def get_number_of_layers(atoms):
    symbols = list(set(atoms.get_chemical_symbols()))
    symbols = sorted(symbols)
    nlayers = []

    for symbol in symbols:
        subatoms = Atoms(list(filter(lambda x: x.symbol == symbol, atoms)))
        pos  = subatoms.positions
        zpos = np.round(pos[:, 2], decimals=4)
        nlayers.append(len(list(set(zpos))))

    return nlayers


def set_tags_by_z(atoms):
    import pandas as pd

    pbc  = atoms.get_pbc()
    cell = atoms.get_cell()

    newatoms = Atoms()
    symbols = list(set(atoms.get_chemical_symbols()))
    symbols = sorted(symbols)

    for symbol in symbols:
        subatoms = Atoms(list(filter(lambda x: x.symbol == symbol, atoms)))
        pos  = subatoms.positions
        zpos = np.round(pos[:, 2], decimals=1)
        bins = list(set(zpos))
        bins = np.sort(bins)
        bins = np.array(bins) + 1.0e-2
        bins = np.insert(bins, 0, 0)

        labels = []
        for i in range(len(bins)-1):
            labels.append(i)

        # tags = pd.cut(zpos, bins=bins, labels=labels).to_list()
        tags = pd.cut(zpos, bins=bins, labels=labels).tolist()

        subatoms.set_tags(tags)
        newatoms += subatoms

    # restore
    newatoms.set_pbc(pbc)
    newatoms.set_cell(cell)

    return newatoms


def remove_layers(atoms=None, element=None, n_layers=1):
    """
    Remove layers of symbol at high-in-z.

    Args:
        atoms (Atoms): Atoms object
        element (str): Element symbol
        n_layers(int): Number of layers (of specified element) to remove
    """
    pbc  = atoms.get_pbc()
    cell = atoms.get_cell()

    atoms_copy = atoms.copy()
    atoms_copy = sort_atoms_by(atoms_copy, xyz="z")  # sort
    atoms_copy = set_tags_by_z(atoms_copy)  # set tags

    newatoms = Atoms()

    tags = atoms_copy.get_tags()
    cond = [i == element for i in atoms_copy.get_chemical_symbols()]

    maxtag = max(list(tags[cond]))

    for i, atom in enumerate(atoms_copy):
        if atom.tag >= maxtag - n_layers + 1 and atom.symbol == element:
            # remove this atom
            pass
        else:
            newatoms += atom

    newatoms.set_pbc(pbc)
    newatoms.set_cell(cell)

    return newatoms


def fix_lower_surface(atoms, adjust_layer=None):
    """
    Fix lower surface atoms. By default, lower half (controled by tag) is fixed.

    Args:
        atoms (Atoms): Atoms object
        adjust_layer (list): List of element-wise layers to adjust fixing. Positive means more layers are fixed.
    """
    from ase.constraints import FixAtoms

    newatoms = atoms.copy()
    newatoms = sort_atoms_by(newatoms, xyz="z")  # sort
    newatoms = set_tags_by_z(newatoms)  # set tags

    # prepare symbol dict
    symbols_ = list(set(atoms.get_chemical_symbols()))
    symbols_ = sorted(symbols_)
    symbols = {}
    for i, sym in enumerate(symbols_):
        symbols.update({sym: i})

    # Determine fixlayer, which is a list of elements. Half of nlayers.
    nlayers = get_number_of_layers(newatoms)

    # check
    div = [i // 2 for i in nlayers]
    mod = [i % 2 for i in nlayers]

    fixlayers = [i + j for (i, j) in zip(div, mod)]

    if adjust_layer is not None:
        fixlayers = [sum(x) for x in zip(fixlayers, adjust_layer)]

    fixlist = []  # list of fixed atoms

    # tags = newatoms.get_tags()
    # minind = np.argmin(tags)
    # maxind = np.argmax(tags)

    # lowest_z  = newatoms[minind].position[2]
    # highest_z = newatoms[maxind].position[2]
    # z_thre = (highest_z - lowest_z) / 2 + lowest_z

    for iatom in newatoms:
        ind = symbols[iatom.symbol]
        # z_pos = iatom.position[2]

        # if iatom.tag < fixlayers[ind] and z_pos < z_thre:
        if iatom.tag < fixlayers[ind]:
            fixlist.append(iatom.index)
        else:
            pass

    constraint = FixAtoms(indices=fixlist)
    newatoms.constraints = constraint

    return newatoms


def find_highest(json, score):
    import pandas as pd

    df = pd.read_json(json)
    df = df.set_index("unique_id")
    df = df.dropna(subset=[score])
    df = df.sort_values(score, ascending=False)

    best = df.iloc[0].name

    return best


def make_step(atoms):
    newatoms = atoms.copy()
    newatoms = sort_atoms_by(newatoms, xyz="z")

    nlayer = get_number_of_layers(newatoms)
    nlayer = nlayer[0]
    perlayer  = len(newatoms) // nlayer
    toplayer  = newatoms[-perlayer:]
    top_layer = sort_atoms_by(toplayer, xyz="y")

    # first remove top layer then add sorted top layer
    del newatoms[-perlayer:]
    newatoms += top_layer

    remove = perlayer // 2

    nstart = perlayer*(nlayer-1)  # index for the atom starting the top layer
    del newatoms[nstart:nstart+remove]

    return newatoms


def mirror_invert(atoms, direction="x"):
    """
    Mirror invert the surface in the specified direction.

    Args:
        atoms: Atoms object
        direction: "x", "y", or "z"
    """
    pos  = atoms.get_positions()
    cell = atoms.cell

    # set position and cell
    if direction == "x":
        pos[:, 0] = -pos[:, 0]
        cell = [[-cell[i][0], cell[i][1], cell[i][2]] for i in range(3)]
    elif direction == "y":
        pos[:, 1] = -pos[:, 1]
        cell = [[cell[i][0], -cell[i][1], cell[i][2]] for i in range(3)]
    elif direction == "z":
        highest_z = pos[:, 2].max()
        atoms.translate([0, 0, -highest_z])
        pos[:, 2] = -pos[:, 2]
        cell = [[cell[i][0], cell[i][1], -cell[i][2]] for i in range(3)]
    else:
        print("direction must be x, y, or z")
        quit()

    atoms.set_positions(pos)

    cell = np.array(cell)
    cell = np.round(cell + 1.0e-5, decimals=4)
    atoms.set_cell(cell)

    return atoms
