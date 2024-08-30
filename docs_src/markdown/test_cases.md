# Test cases
## 1. Oxygen reduction reaction (ORR)
* The ORR is a key reaction in the cathode of H2-O2 fuel cells.
* The overall reaction for the ORR is:

  $$
    \begin{align*}
      & \ce{O2 + 4H+ + 4e-  -> 2H2O}  \ \text{(acidic)} \\
      & \ce{O2 + 2H2O + 4e- -> 4OH-} \ \text{(basic)}
    \end{align*}
  $$

* We will consider the acidic ORR in this test case.
* The acidic ORR consists of the following four elementary reacions;

  $$
    \begin{align*}
      & 1)\ \ce{O2(g) + $*$ + H+ + e- -> OOH$*$} \\
      & 2)\ \ce{OOH$*$ + H+ + e- -> O$*$ + H2O} \\
      & 3)\ \ce{O$*$ + H+ + e- -> OH$*$} \\
      & 4)\ \ce{OH$*$ + H+ + e- -> H2O + $*$}
    \end{align*}
  $$

* Here we will evaluate the Gibbs free energy change ($\Delta G$)
  for each elemeantary reaction, and then take the maximum of these $\Delta G$s.
  This maximum $\Delta G$ is the overpotential ($\eta$) for the ORR.

  $$
    \eta = \max[\Delta G_1, \Delta G_2, \Delta G_3, \Delta G_4]
  $$

* eta is a key parameter to evaluate the catalytic activity, since it is
  the potential difference from the thermodynamically ideal potential (1.23 V).
* To peform the above procedure, we need to evaluate the deltaGs
  for each elementary reactions. This is done by `calc_reaction_energy.`
* DeltaGs should be passed to `rate_oer_and_orr` then $\eta$ is returned.

```python
from microkinetics_toolkit import calc_reaction_energy
from microkinetics_toolkit import rate_oer_and_orr

reaction_file = "orr.txt"
deltaGs = calc_reaction_energy(reaction_file=reaction_file)
eta = rate_oer_and_orr(deltaGs)

print(f"overpotential is {eta} V")
```
