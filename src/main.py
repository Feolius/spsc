import spsc_core
import spsc_data_generator
import scipy.special as sp
import numpy as np
import matplotlib.pyplot as plt

# state = spsc_core.StateSimple.import_file("data/single_well26nm.yml")
# solver = spsc_core.ShrodSolverSimple(state)
# solver.solve()

state = spsc_data_generator.single_well_sloped(26, 10 ** 16)
pot = state.electron_states[0].static_potential
pot.instant_plot()

