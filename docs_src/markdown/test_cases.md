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

* $\eta$ is a key parameter to evaluate the catalytic activity, since it is
  the potential difference from the thermodynamically ideal potential (1.23 V).
* To peform the above procedure, we need to evaluate the $\Delta Gs$.
  for each elementary reactions. This is done by `get_reaction_energy.`
* $\Delta Gs$ should be passed to `get_overpotential_oer_orr` then $\eta$ is returned.

```python
import numpy as np
from microkinetics_toolkit.utils import make_surface_from_cif
from microkinetics_toolkit.utils import remove_layers
from microkinetics_toolkit.utils import replace_element
from microkinetics_toolkit.utils import fix_lower_surface
from microkinetics_toolkit.get_reaction_energy import get_reaction_energy
from microkinetics_toolkit.orr_and_oer import get_overpotential_oer_orr 

# read cif file and make surface
cif_file = "LaMnO3.cif"
surface = make_surface_from_cif(cif_file, indices=(0, 0, 1), vacuum=10.0)

# adjust surface structure
surface = remove_layers(surface, element="La", n_layers=1)
surface = remove_layers(surface, element="O", n_layers=2)

# fix lower part
surface = fix_lower_surface(surface)

# reaction_file = "orr_alkaline.txt"
reaction_file = "orr_alkaline2.txt"

# do calculation
deltaEs = get_reaction_energy(reaction_file=reaction_file, surface=surface, calculator="vasp")
eta = get_overpotential_oer_orr(reaction_file=reaction_file, deltaEs=deltaEs, reaction_type="orr")

print(f"overpotential = {eta:5.3f} eV")
```
