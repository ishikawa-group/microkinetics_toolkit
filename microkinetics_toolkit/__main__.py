def calc_reaction_energy():
    import argparse
    import os
    import numpy as np
    from ase import Atom, Atoms
    from ase.calculators.emt import EMT
    from ase.db import connect
    from ase.io import read, write
    from ase.visualize import view
    from ase.build import surface, add_adsorbate
    from microkinetics_toolkit.utils import get_number_of_reaction
    from microkinetics_toolkit.utils import get_reac_and_prod
    from microkinetics_toolkit.utils import get_adsorbate_type
    from microkinetics_toolkit.utils import make_surface_from_cif

    parser = argparse.ArgumentParser()
    parser.add_argument("--reaction_file", default="oer.txt")
    parser.add_argument("--cif_file", default="pdo.cif")

    args = parser.parse_args()
    reaction_file = args.reaction_file
    cif_file = args.cif_file

    r_ads, r_site, r_coef, p_ads, p_site, p_coef = get_reac_and_prod(reaction_file)
    rxn_num = get_number_of_reaction(reaction_file)

    # load molecule collection from ase
    database = connect("data/g2plus.json")

    # reaction energy
    deltaEs = np.array([])

    # define calculator for molecules and surfaces separately
    calc_mol = EMT()
    calc_surf = EMT()

    for irxn in range(rxn_num):
        print(f"irxn = {irxn}")

        energies = {"reactant": 0.0, "product": 0.0}

        for side in ["reactant", "product"]:
            if side == "reactant":
                mols, sites, coefs = r_ads[irxn], r_site[irxn], r_coef[irxn]
            elif side == "product":
                mols, sites, coefs = p_ads[irxn], p_site[irxn], p_coef[irxn]
            else:
                print("some error")
                quit()

            E = 0.0

            surf = make_surface_from_cif(cif_file)

            for imol, mol in enumerate(mols):
                if mol[0] == "surf":
                    atoms = surf
                else:
                    id_ = database.get(name=mol[0]).id
                    atoms = database.get_atoms(id=id_)

                site = sites[imol][0]
                ads_type = get_adsorbate_type(atoms, site)

                rotation = {"HO" : [180, "x"],
                            "HO2": [180, "x"]}

                if ads_type == "gaseous":
                    atoms.calc = calc_mol
                elif ads_type == "surface":
                    atoms.calc = calc_surf
                elif ads_type == "adsorbed":
                    adsorbate = atoms
                    formula = adsorbate.get_chemical_formula()
                    if formula in rotation:
                        adsorbate.rotate(*rotation[formula])

                    height = 1.8
                    position = adsorbate.positions[0][:2]
                    add_adsorbate(surf, adsorbate, offset=(0, 0), position=position, height=height)

                    atoms = surf.copy()
                    atoms.calc = calc_surf
                    # view(atoms)
                else:
                    print("some error")
                    quit()

                energy = atoms.get_potential_energy()
                E += coefs[imol]*energy

            energies[side] = E

        deltaE  = energies["product"] - energies["reactant"]
        deltaEs = np.append(deltaEs, deltaE)

    return deltaEs


def calc_rate(deltaEs):
    rate = 0
    return rate


if __name__ == "__main__":
    deltaEs = calc_reaction_energy()
    print(f"deltaEs = {deltaEs}")
    rate = calc_rate(deltaEs)
    print(f"rate = {rate}")
