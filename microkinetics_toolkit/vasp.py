def set_vasp_calculator(atom_type="molecule", dfttype="gga", kpt=1, do_optimization=False):
    """
    Set up a calculator using VASP.

    Args:
        atom_type: "molecule", "surface" or "solid".
        kpt: k-points in x and y directions.
        do_optimization: True or False.
    """
    import json
    from ase.calculators.vasp import Vasp

    if atom_type == "molecule":
        kpts = [1, 1, 1]
        ismear = 0
        lreal = True
        ldipol = False
        idipol = None
    elif atom_type == "surface":
        kpts = [kpt, kpt, 1]
        ismear = 0
        lreal = True
        ldipol = True
        idipol = 3
    elif atom_type == "solid":
        kpts = [kpt, kpt, kpt]
        ismear = 0
        lreal = False
        ldipol = False
        idipol = None
    else:
        print("some error")
        quit()

    # common setting
    xc = "pbe"
    encut = 500.0
    ediff  = 1.0e-4
    ediffg = -10.0e-2
    lorbit = 10
    algo = "Normal"
    nelmin = 5
    nelm = 30
    npar = 4
    nsim = npar
    ispin = 2
    kgamma = True
    setups = {"Cr": "_pv", "Mn": "_pv", "Fe": "_pv"}
    lasph = True

    # DFT + U
    if dfttype == "plus_u":
        ldau = True
        ldautype = 2
        u_param_file = "data/u_parameter.json"
        with open(u_param_file) as f:
            ldau_luj = json.load(f)
    else:
        ldau = None
        ldautype = None
        ldau_luj = None

    # meta-gga
    if dfttype == "meta-gga":
        metagga = "r2scan"
    else:
        metagga = None

    # geometry optimization related
    if do_optimization:
        ibrion = 2
        potim = 0.1
        nsw = 1
    else:
        ibrion = 0
        potim = 0.0
        nsw = 0

    calc = Vasp(prec="Normal", xc=xc, metagga=metagga, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, ldipol=ldipol, idipol=idipol, setups=setups, lasph=True,
                ispin=ispin, npar=npar, nsim=nsim, nelmin=nelmin, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj,
                )

    return calc


def set_calc_directory(atoms=None, dirname=None, formula=None):
    directory = "work_" + dirname + "/" + formula
    atoms.calc.directory = directory
    return None


def set_lmaxmix(atoms=None):
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

    return None
