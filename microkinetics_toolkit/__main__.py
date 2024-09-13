if __name__ == "__main__":
    import numpy as np
    from microkinetics_toolkit.utils import make_surface_from_cif
    from microkinetics_toolkit.utils import remove_layers
    from microkinetics_toolkit.utils import replace_element
    from microkinetics_toolkit.utils import fix_lower_surface
    from microkinetics_toolkit.calc_reaction_energy import calc_reaction_energy
    from microkinetics_toolkit.orr_and_oer import calc_overpotential_oer_orr 
    from ase.visualize import view

    cif_file = "LaMnO3.cif"
    surface = make_surface_from_cif(cif_file, indices=(0, 0, 1), vacuum=10.0)

    # for EMT
    use_emt = False
    if use_emt:
        surface = replace_element(surface, from_element="La", to_element="Al")
        surface = replace_element(surface, from_element="Mn", to_element="Pt")
        surface = remove_layers(surface, element="Al", n_layers=1)
    else:
        surface = remove_layers(surface, element="La", n_layers=1)

    surface = remove_layers(surface, element="O", n_layers=2)
    surface = fix_lower_surface(surface)

    # reaction_file = "orr2.txt"
    # reaction_file = "orr_alkaline.txt"
    reaction_file = "orr_alkaline2.txt"

    deltaEs = calc_reaction_energy(reaction_file=reaction_file, surface=surface, calculator="vasp", verbose=True)
    eta = calc_overpotential_oer_orr(reaction_file=reaction_file, deltaEs=deltaEs, reaction_type="orr", verbose=True)
    eta = np.abs(eta)

    print(f"eta = {eta:5.3f} eV")
