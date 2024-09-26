def get_overpotential_oer_orr(reaction_file, deltaEs, T=298.15, reaction_type="oer", energy_shift=None, verbose=False):
    """
    Calculate overpotential for OER or ORR.
    """
    import numpy as np
    import matplotlib.pyplot as plt
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

    if energy_shift is not None:
        deltaGs += np.array(energy_shift)

    if verbose:
        print(f"max of deltaGs = {np.max(deltaGs):5.3f} eV")

    if reaction_type == "oer":
        eta = np.max(deltaGs) - 1.23
    elif reaction_type == "orr":
        eta = 1.23 - np.max(deltaGs)
        eta = np.abs(eta)  # necessary?
    else:
        print("some error")
        quit()

    np.set_printoptions(precision=3)

    print(f"deltaGs = {deltaGs}")

    # make deltaG relative
    deltaGs_rel = deltaGs - deltaGs[0]
    deltaGs_rel = np.append(deltaGs_rel, 0)

    print(f"deltaGs_rel = {deltaGs_rel}")

    # plot
    fig_name = "test.png"
    plt.plot(deltaGs_rel, "o")
    plt.savefig(fig_name)

    return eta
