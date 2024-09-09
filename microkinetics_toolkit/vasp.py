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
        lreal = "True"
    elif atom_type == "surface":
        kpts = [kpt, kpt, 1]
        ismear = 0 if kpt == 1 else 1
        lreal = "False"
    elif atom_type == "solid":
        kpts = [kpt, kpt, kpt]
        ismear = 0 if kpt == 1 else 1
        lreal = "False"
    else:
        print("some error")
        quit()

    # common setting
    xc = "pbe"
    encut = 400.0
    ediff  =  1.0e-3
    ediffg = -5.0e-1
    lorbit = 10
    algo = "Normal"
    nelm = 30
    npar = 2
    nsim = npar
    ispin = 2
    kgamma = True

    # DFT + U
    if dfttype == "plus_u":
        ldau = True
        ldautype = 2

        u_param_file = "data/u_parameter.json"
        with open(u_param_file) as f:
            ldau_luj = json.load(f)

    # geometry optimization related
    if do_optimization:
        ibrion = 2
        potim = 0.1
        nsw = 100
    else:
        ibrion = 0
        potim = 0.0
        nsw = 0

    calc = Vasp(prec="Normal", xc=xc, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, 
                ispin=ispin, npar=npar, nsim=nsim, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                )

    if dfttype == "plus_u":
        calc = Vasp(prec="Normal", xc=xc, encut=encut, kpts=kpts, ismear=ismear, ediff=ediff, ediffg=ediffg,
                    ibrion=ibrion, potim=potim, nsw=nsw, algo=algo, 
                    ispin=ispin, npar=npar, nsim=nsim, nelm=nelm, lreal=lreal, lorbit=lorbit, kgamma=kgamma,
                    ldau=ldau, ldautype=ldautype, ldau_luj=ldau_luj,
                    )

    return calc

