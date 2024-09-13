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
<<<<<<< HEAD
<<<<<<< HEAD
        lreal = True
        ldipol = False
        idipol = None
    elif atom_type == "surface":
        kpts = [kpt, kpt, 1]
        ismear = 0
        lreal = True  # True or False
        ldipol = True
        idipol = 3
    elif atom_type == "solid":
        kpts = [kpt, kpt, kpt]
        ismear = 0
        lreal = False
        ldipol = False
        idipol = None
=======
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
        lreal = "True"
    elif atom_type == "surface":
        kpts = [kpt, kpt, 1]
        ismear = 0 if kpt == 1 else 1
        lreal = "False"
    elif atom_type == "solid":
        kpts = [kpt, kpt, kpt]
        ismear = 0 if kpt == 1 else 1
        lreal = "False"
<<<<<<< HEAD
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    else:
        print("some error")
        quit()

    # common setting
    xc = "pbe"
<<<<<<< HEAD
<<<<<<< HEAD
    encut = 500.0
=======
    encut = 400.0
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
    encut = 400.0
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    ediff  =  1.0e-4
    ediffg = -1.0e-0
    lorbit = 10
    algo = "Normal"
    nelmin = 5
    nelm = 30
<<<<<<< HEAD
<<<<<<< HEAD
    npar = 4
    nsim = npar
    ispin = 2
    kgamma = True
    setups = {"Mn": "_pv", "Fe": "_pv"}
    lasph = True
=======
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    npar = 2
    nsim = npar
    ispin = 2
    kgamma = True
<<<<<<< HEAD
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878

    # DFT + U
    if dfttype == "plus_u":
        ldau = True
        ldautype = 2
<<<<<<< HEAD
<<<<<<< HEAD
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
=======
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878

        u_param_file = "data/u_parameter.json"
        with open(u_param_file) as f:
            ldau_luj = json.load(f)
<<<<<<< HEAD
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878

    # geometry optimization related
    if do_optimization:
        ibrion = 2
        potim = 0.1
<<<<<<< HEAD
<<<<<<< HEAD
        nsw = 10
=======
        nsw = 100
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
        nsw = 100
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    else:
        ibrion = 0
        potim = 0.0
        nsw = 0

<<<<<<< HEAD
<<<<<<< HEAD
    calc = Vasp(prec="Normal", xc=xc, metagga=metagga, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, ldipol=ldipol, idipol=idipol, setups=setups, lasph=True,
                ispin=ispin, npar=npar, nsim=nsim, nelmin=nelmin, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj,
                )

=======
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    calc = Vasp(prec="Normal", xc=xc, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, 
                ispin=ispin, npar=npar, nsim=nsim, nelmin=nelmin, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                )

    if dfttype == "plus_u":
        calc = Vasp(prec="Normal", xc=xc, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                    ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, 
                    ispin=ispin, npar=npar, nsim=nsim, nelmin=nelmin, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                    ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj,
                    )

<<<<<<< HEAD
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
=======
>>>>>>> d14e415f7deb6fdc8dd8c6128555cf0d4dd05878
    return calc

