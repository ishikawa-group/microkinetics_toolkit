def get_reaction_energy(reaction_file="oer.txt", surface=None, calculator="emt", verbose=False):
    """
    Calculate reaction energy for each reaction.
    """
    import os
    import numpy as np
    from ase.build import add_adsorbate
    from ase.calculators.emt import EMT
    from ase.db import connect
    from ase.visualize import view
    from microkinetics_toolkit.utils import get_adsorbate_type
    from microkinetics_toolkit.utils import get_number_of_reaction
    from microkinetics_toolkit.utils import get_reac_and_prod
    from microkinetics_toolkit.vasp import set_vasp_calculator

    r_ads, r_site, r_coef, p_ads, p_site, p_coef = get_reac_and_prod(reaction_file)
    rxn_num = get_number_of_reaction(reaction_file)

    # load molecule collection
    database_file = "data/g2plus.json"
    if os.path.exists(database_file):
        database = connect(database_file)
    else:
        raise FileNotFoundError

    # temporary database
    tmpdbfile = "tmp.db"
    tmpdb = connect(tmpdbfile)

    # reaction energy
    deltaEs = np.array([])

    # define calculator for molecules and surfaces separately
    if "emt" in calculator:
        calc_mol  = EMT()
        calc_surf = EMT()
    elif "vasp" in calculator:
        calc_mol  = set_vasp_calculator(atom_type="molecule", do_optimization=True, dfttype="plus_u", kpt=2)
        calc_surf = set_vasp_calculator(atom_type="surface", do_optimization=True, dfttype="plus_u", kpt=2)
    elif "ocp" in valculator:
        calc_mol  = set_ocp_calculator()  # do not work
        calc_surf = set_ocp_calculator()
    else:
        raise ValueError("Choose from emt, vasp, ocp.")

    # rotational angle for adsorbed molecules
    rotation = {"HO": [180, "x"], "HO2": [180, "x"], "O2": [90, "x"]}

    for irxn in range(rxn_num):
        energies = {"reactant": 0.0, "product": 0.0}

        for side in ["reactant", "product"]:
            surface_ = surface.copy()

            # assume high spin for surface
            symbols = surface_.get_chemical_symbols()
            init_magmom = [1.0 if x == "Mn" else 0.0 for x in symbols]
            surface_.set_initial_magnetic_moments(init_magmom)

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
                    try:
                        id_ = database.get(name=mol[0]).id
                    except KeyError:
                        print(f"{mol[0]} not found in {database_file}")
                        quit()
                    else:
                        atoms = database.get_atoms(id=id_)

                site = sites[imol][0]
                ads_type = get_adsorbate_type(atoms, site)

                if ads_type == "gaseous":
                    if mol[0] == "surf":
                        atoms.calc = calc_surf
                    else:
                        atoms.calc = calc_mol
                        atoms.cell = [20, 20, 20]
                        atoms.center()
                elif ads_type == "surface":
                    atoms.calc = calc_surf
                elif ads_type == "adsorbed":
                    adsorbate = atoms
                    tmp = adsorbate.get_chemical_formula()

                    if tmp in rotation:
                        adsorbate.rotate(*rotation[tmp])

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

                # Identification done. Look for temporary database for identical system.
                formula = atoms.get_chemical_formula()

                try:
                    previous = tmpdb.get(name=formula)
                except KeyError:
                    first_time = True
                    if verbose:
                        print(f"Calculating {formula} ... new calculation", flush=True)
                else:
                    first_time = False
                    if verbose:
                        print(f"Calculating {formula} ... reuse previous value", flush=True)

                if not first_time:
                    #
                    # use previous energy
                    #
                    energy = tmpdb.get(id=previous.id).data.energy
                else:
                    #
                    # perform energy calculation
                    #

                    # final setting before calcuation
                    atoms.pbc = True
                    directory = "work" + "/" + formula
                    atoms.calc.directory = directory

                    # lmaxmix setting
                    symbols = atoms.get_chemical_symbols()
                    d_elements = ["Mo", "Fe", "Cr"]
                    f_elements = ["La", "Ce", "Pr", "Nd"]

                    if len(set(symbols) & set(d_elements)) != 0:
                        # has some d_elements
                        atoms.calc.set(lmaxmix=4)

                    if len(set(symbols) & set(f_elements)) != 0:
                        # has some f_elements
                        atoms.calc.set(lmaxmix=6)

                    # do calculation
                    energy = atoms.get_potential_energy()

                E += coefs[imol]*energy

                # recording to database
                if first_time:
                    id = tmpdb.reserve(name=formula)
                    if id is None:
                        # somebody is writing to db
                        continue
                    else:
                        # ready to write
                        name = formula
                        tmpdb.write(atoms, name=name, id=id, data={"energy": energy})

            energies[side] = E

        deltaE  = energies["product"] - energies["reactant"]
        deltaEs = np.append(deltaEs, deltaE)

    return deltaEs

