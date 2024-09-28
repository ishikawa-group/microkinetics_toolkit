def register(db=None, atoms=None, formula=None, data=None):
    if db is None:
        print("in registar: no database is given")
        quit()

    if atoms is None:
        print("in registar: nothing to registar")
        quit()

    id = db.reserve(name=formula)
    if id is not None:
        name = formula
        db.write(atoms, name=name, id=id, data=data)


def get_past_atoms(db=None, atoms=None):
    if db is None:
        print("in get_past_data: no database is given")
        quit()

    if atoms is None:
        print("in get_past_data: nothing to registar")
        quit()

    formula = atoms.get_chemical_formula()
    try:
        past = tmpdb.get(name=surf_formula)
    except:
        # not found - return atoms as is
        first = True
        return atoms, first
    else:
        # found - return old atoms
        first = False
        atoms_ = tmpdb.get_atoms(id=past.id).copy()
        return atoms_, first


def get_reaction_energy(reaction_file="oer.txt", surface=None, calculator="emt", verbose=False, dirname=None):
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
    from microkinetics_toolkit.vasp import set_directory
    from microkinetics_toolkit.vasp import set_lmaxmix

    r_ads, r_site, r_coef, p_ads, p_site, p_coef = get_reac_and_prod(reaction_file)
    rxn_num = get_number_of_reaction(reaction_file)

    # load molecule collection
    database_file = "data/g2plus.json"
    if os.path.exists(database_file):
        database = connect(database_file)
    else:
        raise FileNotFoundError

    # temporary database
    tmpdbfile = "tmp_" + dirname + ".db"
    tmpdb = connect(tmpdbfile)

    # reaction energy
    deltaEs = np.array([])

    # define calculator for molecules and surfaces separately
    if "emt" in calculator:
        calc_mol  = EMT()
        calc_surf = EMT()
    elif "vasp" in calculator:
        kpt = 1
        dfttype = "gga"  # "plus_u"
        calc_mol  = set_vasp_calculator(atom_type="molecule", do_optimization=True, dfttype=dfttype, kpt=kpt)
        calc_surf = set_vasp_calculator(atom_type="surface", do_optimization=True, dfttype=dfttype, kpt=kpt)
    elif "ocp" in valculator:
        calc_mol  = set_ocp_calculator()  # do not work
        calc_surf = set_ocp_calculator()
    else:
        raise ValueError("Choose from emt, vasp, ocp.")

    # rotational angle for adsorbed molecules
    # rotation = {"HO": [180, "x"], "HO2": [180, "x"], "O2": [90, "x"]}
    rotation = {"HO": [160, "x"], "HO2": [160, "x"], "O2": [70, "x"]}

    # spin-polarized or not for adsorbed molecules
    closed_shell_molecules = ["H2", "HO", "H2O"]

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
                        if mol[0] in closed_shell_molecules:
                            atoms.calc.set(ispin=1)

                elif ads_type == "surface":
                    atoms.calc = calc_surf

                elif ads_type == "adsorbed":
                    adsorbate = atoms
                    tmp = adsorbate.get_chemical_formula()

                    if tmp in rotation:
                        adsorbate.rotate(*rotation[tmp])

                    height = 1.8
                    offset = (0.0, 0.25)  # (0.0, 0.50) for smallest cell
                    position = adsorbate.positions[0][:2]

                    # get surf part from tmpdb of ads + surf case, as it should be done beforehand
                    surf_formula = surface_.get_chemical_formula()
                    surface_, first = get_past_atoms(db=tmpdb, atoms=surface_)
                    if first:
                        print("first time for bare surface - do calculation here")
                        directory = "work_" + dirname + "/" + formula
                        set_calc_directory(atoms=atoms, directory=directory)
                        set_lmaxmix(atoms=atoms)
                        energy = atoms.get_potential_energy()
                        register(db=tmpdb, atoms=atoms, formula=formula, data={"energy": energy})
                    else:
                        print("found data for bare surface")
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
                    # use previous energy
                    energy = tmpdb.get(id=previous.id).data.energy
                else:
                    # perform energy calculation
                    directory = "work_" + dirname + "/" + formula
                    set_calc_directory(atoms=atoms, directory=directory)
                    set_lmaxmix(atoms=atoms)
                    energy = atoms.get_potential_energy()

                E += coefs[imol]*energy

                # recording to database
                if first_time:
                    register(db=tmpdb, atoms=atoms, formula=formula, data={"energy": energy})

                """
                if first_time:
                    id = tmpdb.reserve(name=formula)
                    if id is None:
                        # somebody is writing to db
                        continue
                    else:
                        # ready to write
                        name = formula
                        tmpdb.write(atoms, name=name, id=id, data={"energy": energy})
                """

            energies[side] = E

        deltaE  = energies["product"] - energies["reactant"]
        deltaEs = np.append(deltaEs, deltaE)

    return deltaEs
