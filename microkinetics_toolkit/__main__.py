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


def calc_reaction_energy(reaction_file="oer.txt", cif_file="pdo.cif"):
    """
    Calculate reaction energy for each reaction.
    """
    r_ads, r_site, r_coef, p_ads, p_site, p_coef = get_reac_and_prod(reaction_file)
    rxn_num = get_number_of_reaction(reaction_file)

    # load molecule collection from ase
    database = connect("data/g2plus.json")

    # reaction energy
    deltaEs = np.array([])

    # define calculator for molecules and surfaces separately
    calc_mol = EMT()
    calc_surf = EMT()

    # calc_mol  = set_ocp_calculator()  # do not work
    # calc_surf = set_ocp_calculator()

    # rotationa angle for adsorbed molecules
    rotation = {"HO": [180, "x"], "HO2": [180, "x"]}

    for irxn in range(rxn_num):
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


def calc_overpotential_oer_orr(reaction_file, deltaEs, T=298.15):
    """
    Calculate overpotential for OER or ORR.
    """
    rxn_num = get_number_of_reaction(reaction_file)

    zpe = {"H2": 0.0, "H2O": 0.0, "OHads": 0.0, "Oads": 0.0, "OOHads": 0.0}
    S = {"H2": 0.0, "H2O": 0.0, "O2": 0.0}

    # ZPE in eV
    zpe["H2"]  = 0.27
    zpe["H2O"] = 0.56
    zpe["OHads"] = 0.36
    zpe["Oads"] = 0.07
    zpe["OOHads"] = 0.40
    zpe["O2"] = 0.05*2

    # entropy in eV/K
    S["H2"] = 0.41/T
    S["H2O"] = 0.67/T
    S["O2"] = 0.32*2/T

    # loss in entropy in each reaction
    deltaSs = np.zeros(rxn_num)
    deltaSs[0] = 0.5*S["H2"] - S["H2O"]
    deltaSs[1] = 0.5*S["H2"]
    deltaSs[2] = 0.5*S["H2"] - S["H2O"]
    deltaSs[3] = 2.0*S["H2O"] - 1.5*S["H2"]

    deltaZPEs = np.zeros(rxn_num)
    deltaZPEs[0] = zpe["OHads"] + 0.5*zpe["H2"] - zpe["H2O"]
    deltaZPEs[1] = zpe["Oads"] + 0.5*zpe["H2"] - zpe["OHads"]
    deltaZPEs[2] = zpe["OOHads"] + 0.5*zpe["H2"] - zpe["Oads"] - zpe["H2O"]
    deltaZPEs[3] = 2.0*zpe["H2O"] - 1.5*zpe["H2"] - zpe["OOHads"]

    deltaEs = np.array(deltaEs)
    deltaHs = deltaEs + deltaZPEs
    deltaGs = deltaHs - T*deltaSs

    eta_oer = np.max(deltaGs) - 1.23

    return eta_oer


def set_vasp_calculator():
    """
    Set up a calculator using VASP.
    """
    return calc


def set_ocp_calculator():
    """
    Set up a calculator using Neural Network Potential with OCP.
    """
    from fairchem.core.models.model_registry import available_pretrained_models
    from fairchem.core.models.model_registry import model_name_to_local_file
    from fairchem.core.common.relaxation.ase_utils import OCPCalculator

    checkpoint_path = model_name_to_local_file("GemNet-dT-S2EFS-OC22", local_cache="./downloaded_checkpoints/")
    calc = OCPCalculator(checkpoint_path=checkpoint_path)

    return calc


if __name__ == "__main__":
    from ase.build import fcc111

    deltaEs = calc_reaction_energy(reaction_file="oer.txt", cif_file="pdo.cif")
    print(f"deltaEs = {deltaEs}")
    eta = calc_overpotential_oer_orr(reaction_file="oer.txt", deltaEs=deltaEs)
    print(f"eta = {eta:5.3f} eV")
