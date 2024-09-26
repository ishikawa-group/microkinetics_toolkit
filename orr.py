if __name__ == "__main__":
    import numpy as np
    import argparse
    from microkinetics_toolkit.utils import make_surface_from_cif
    from microkinetics_toolkit.utils import remove_layers
    from microkinetics_toolkit.utils import replace_element
    from microkinetics_toolkit.utils import fix_lower_surface
    from microkinetics_toolkit.get_reaction_energy import get_reaction_energy
    from microkinetics_toolkit.orr_and_oer import get_overpotential_oer_orr 
    from ase.visualize import view

    parser = argparse.ArgumentParser()
    parser.add_argument("--unique_id", default="0")
    parser.add_argument("--replace_percent", default=0)
    args = parser.parse_args()
    unique_id = args.unique_id
    replace_percent = int(args.replace_percent)

    cif_file = "LaMnO3.cif"
    surface = make_surface_from_cif(cif_file, indices=[0, 0, 1], repeat=[2, 2, 1], vacuum=7.0)

    # make random replacement
    surface = replace_element(surface, from_element="Mn", to_element="Cr", percent=100)
    # surface = replace_element(surface, from_element="Fe", to_element="Cr", percent=replace_percent)

    surface = remove_layers(surface, element="La", n_layers=1)
    surface = remove_layers(surface, element="O", n_layers=2)

    surface = fix_lower_surface(surface)

    # reaction_file = "orr_alkaline.txt"
    reaction_file = "orr_alkaline2.txt"
    # reaction_file = "orr_alkaline3.txt"
    # energy_shift = [-4.92, 0, 0, 0]

    deltaEs = get_reaction_energy(reaction_file=reaction_file, surface=surface, calculator="vasp", verbose=True, dirname=unique_id)
    # eta = get_overpotential_oer_orr(reaction_file=reaction_file, deltaEs=deltaEs, reaction_type="orr", verbose=True, energy_shift=energy_shift)
    eta = get_overpotential_oer_orr(reaction_file=reaction_file, deltaEs=deltaEs, reaction_type="orr", verbose=True)
    eta = np.abs(eta)

    print(f"eta = {eta:5.3f} eV")
