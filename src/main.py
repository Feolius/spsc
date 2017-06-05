import scipy.special as sp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import spsc_data_generator
import spsc_data
import spsc_constants
import spsc_core

# state = spsc_data_generator.single_well_sloped(spsc_data.LengthValue(26, "nm"), spsc_data.ElectricFieldValue(7 * 10 ** -4, "V_per_m"))
# pot = state.electron_states[0].static_potential
# pot.instant_plot()
# solver = spsc_core.ShrodSolverSimple(state)
# solver.solve()
# state = spsc_data_generator.superlattice_well(8, spsc_data.LengthValue(26, "nm"), spsc_data.LengthValue(2.3, "nm"), spsc_data.LengthValue(1.1, "nm"))
state = spsc_data_generator.superlattice_well_sloped(8, spsc_data.LengthValue(26, "nm"),
                                                     spsc_data.LengthValue(2.3, "nm"), spsc_data.LengthValue(1.1, "nm"),
                                                     spsc_data.ElectricFieldValue(10 ** -4, "V_per_m"))
# pot = state.electron_states[0].static_potential
# pot.instant_plot()
# solver = spsc_core.LatticeSymmetrySolver(state)
# solver.solve()




