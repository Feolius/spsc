import spsc_core
import spsc_data_generator
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

# state = spsc_core.StateSimple.import_file("data/single_well26nm.yml")
# solver = spsc_core.ShrodSolverSimple(state)
# solver.solve()

# state = spsc_data_generator.single_well_sloped(26, 1)
# pot = state.electron_states[0].static_potential
# pot.instant_plot()
args = np.arange(-15, 1, 0.1)
airyValues = [sp.airy(x) for x in args]
# airyValues = np.array(airyValues)
airyA = []
airyB = []
for v in airyValues:
    airyA.append(v[0])
    airyB.append(v[2])
plt.gcf().clear()
# plt.plot(args, airyA)
plt.plot(args, airyB)
plt.show()
