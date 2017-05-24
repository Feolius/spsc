import scipy.special as sp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import spsc_data_generator
import spsc_data
import spsc_constants
import spsc_core

# state = spsc_data_generator.single_well(spsc_data.LengthValue(26, "nm"))
# pot = state.electron_states[0].static_potential
# pot.convert_to("erg")
# pot.instant_plot()
# solver = spsc_core.ShrodSolverSimple(state)
# solver.solve()

state = spsc_data_generator.single_well_sloped(spsc_data.LengthValue(26, "nm"), spsc_data.ElectricFieldValue(5 * 10 ** -4, "V_per_m"))
# pot = state.electron_states[0].static_potential
# pot.instant_plot()
solver = spsc_core.ShrodSolverSimple(state)
solver.solve()




