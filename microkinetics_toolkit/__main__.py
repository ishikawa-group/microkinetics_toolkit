def calc_reaction_energy(reaction_file="oer.txt", surface=None, verbose=False):
    """
    Calculate reaction energy for each reaction.
    """
    import numpy as np
    from ase.build import add_adsorbate
    # from ase.calculators.emt import EMT
    from ase.calculators.vasp import VASP
    from ase.db import connect
    from ase.visualize import view
    from microkinetics_toolkit.utils import get_adsorbate_type
    from microkinetics_toolkit.utils import get_number_of_reaction
    from microkinetics_toolkit.utils import get_reac_and_prod

    r_ads, r_site, r_coef, p_ads, p_site, p_coef = get_reac_and_prod(reaction_file)
    rxn_num = get_number_of_reaction(reaction_file)

    # load molecule collection from ase
    database = connect("data/g2plus.json")

    # reaction energy
    deltaEs = np.array([])

    # define calculator for molecules and surfaces separately
    calculator = "emt"
    if "emt" in calculator:
        calc_mol  = EMT()
        calc_surf = EMT()
    elif "vasp" in calculator:
        calc_mol  = set_vasp_calculator()
        calc_surf = set_vasp_calculator()
    ellif "ocp" in valculator:
        calc_mol  = set_ocp_calculator()  # do not work
        calc_surf = set_ocp_calculator()

    else:
        print("emt or vasp")
        quit()


    # rotationa angle for adsorbed molecules
    rotation = {"HO": [180, "x"], "HO2": [180, "x"]}

    for irxn in range(rxn_num):
        energies = {"reactant": 0.0, "product": 0.0}

        for side in ["reactant", "product"]:
            surface_ = surface.copy()

            if side == "reactant":
                mols, sites, coefs = r_ads[irxn], r_site[irxn], r_coef[irxn]
            elif side == "product":
                mols, sites, coefs = p_ads[irxn], p_site[irxn], p_coef[irxn]
            else:
                print("some error")
                quit()

            E = 0.0
            for imol, mol in enumerate(mols):
                if mol[0] == "surf":
                    atoms = surface_
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
                    offset = (0.0, 0.5)
                    position = adsorbate.positions[0][:2]
                    add_adsorbate(surface_, adsorbate, offset=offset, position=position, height=height)

                    atoms = surface_.copy()
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


def calc_overpotential_oer_orr(reaction_file, deltaEs, T=298.15, reaction_type="oer", verbose=False):
    """
    Calculate overpotential for OER or ORR.
    """
    import numpy as np
    from microkinetics_toolkit.utils import get_number_of_reaction

    rxn_num = get_number_of_reaction(reaction_file)

    zpe = {"H2": 0.0, "H2O": 0.0, "OHads": 0.0, "Oads": 0.0, "OOHads": 0.0}
    S = {"H2": 0.0, "H2O": 0.0, "O2": 0.0}

    # ZPE in eV
    zpe["H2"] = 0.27
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
    deltaZPEs = np.zeros(rxn_num)

    reaction_type = reaction_type.lower()
    if reaction_type == "oer":
        deltaSs[0] = 0.5*S["H2"] - S["H2O"]
        deltaSs[1] = 0.5*S["H2"]
        deltaSs[2] = 0.5*S["H2"] - S["H2O"]
        deltaSs[3] = 2.0*S["H2O"] - 1.5*S["H2"]

        deltaZPEs[0] = zpe["OHads"] + 0.5*zpe["H2"] - zpe["H2O"]
        deltaZPEs[1] = zpe["Oads"] + 0.5*zpe["H2"] - zpe["OHads"]
        deltaZPEs[2] = zpe["OOHads"] + 0.5*zpe["H2"] - zpe["Oads"] - zpe["H2O"]
        deltaZPEs[3] = 2.0*zpe["H2O"] - 1.5*zpe["H2"] - zpe["OOHads"]

    elif reaction_type == "orr":
        deltaSs[0] = - S["O2"] - S["H2"]
        deltaSs[1] = S["H2O"] - 0.5*S["H2"]
        deltaSs[2] = - 0.5*S["H2"]
        deltaSs[3] = S["H2O"] - 0.5*S["H2"]

        deltaZPEs[0] = zpe["OOHads"] - 0.5*zpe["H2"] - zpe["O2"]
        deltaZPEs[1] = zpe["Oads"] + zpe["H2O"] - 0.5*zpe["H2"] - zpe["OOHads"]
        deltaZPEs[2] = zpe["OHads"] - 0.5*zpe["H2"] - zpe["Oads"]
        deltaZPEs[3] = zpe["H2O"] - 0.5*zpe["H2"] - zpe["OHads"]

    else:
        print("some error")
        quit()

    deltaEs = np.array(deltaEs)
    deltaHs = deltaEs + deltaZPEs
    deltaGs = deltaHs - T*deltaSs

    if verbose:
        print(f"max of deltaGs = {np.max(deltaGs):5.3f} eV")
    if reaction_type == "oer":
        eta = np.max(deltaGs) - 1.23
    elif reaction_type == "orr":
        eta = 1.23 - np.max(deltaGs)
    else:
        print("some error")
        quit()

    return eta


def set_vasp_calculator():
    """
    Set up a calculator using VASP.
    """
    # return calc
    return None


def set_ocp_calculator():
    """
    Set up a calculator using Neural Network Potential with OCP.
    """
    from fairchem.core.common.relaxation.ase_utils import OCPCalculator
    from fairchem.core.models.model_registry import model_name_to_local_file

    checkpoint_path = model_name_to_local_file("GemNet-dT-S2EFS-OC22", local_cache="./downloaded_checkpoints/")
    calc = OCPCalculator(checkpoint_path=checkpoint_path)

    return calc


if __name__ == "__main__":
    from microkinetics_toolkit.utils import make_surface_from_cif
    from microkinetics_toolkit.utils import remove_layers
    from microkinetics_toolkit.utils import replace_element
    from microkinetics_toolkit.utils import fix_lower_surface
    from ase.visualize import view

    cif_file = "LaMnO3.cif"
    surface = make_surface_from_cif(cif_file, indices=(0, 0, 1), vacuum=10.0)

    # for EMT
    surface = replace_element(surface, from_element="La", to_element="Al")
    surface = replace_element(surface, from_element="Mn", to_element="Pt")

    surface = remove_layers(surface, element="O", n_layers=2)
    surface = remove_layers(surface, element="Al", n_layers=1)
    surface = fix_lower_surface(surface)

    reaction_file = "orr2.txt"
    deltaEs = calc_reaction_energy(reaction_file=reaction_file, surface=surface)
    print(f"deltaEs = {deltaEs}")
    eta = calc_overpotential_oer_orr(reaction_file=reaction_file, deltaEs=deltaEs, reaction_type="orr", verbose=True)
    print(f"eta = {eta:5.3f} eV")
