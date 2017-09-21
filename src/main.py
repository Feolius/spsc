import scipy.special as sp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import spsc_data_generator
import spsc_data
import spsc_constants
import spsc_core
import yaml

# state = spsc_data_generator.single_well(spsc_data.LengthValue(26, "nm"))
# solver = spsc_core.SymmetricWellSolver(state)
# solver.solve()
# state = spsc_data_generator.single_well_sloped(spsc_data.LengthValue(26, "nm"), spsc_data.ElectricFieldValue(7 * 10 ** -4, "V_per_m"))
# solver = spsc_core.SlopedWellSolver(state)
# solver.solve()
state = spsc_data_generator.x_electrons_superlattice('data/data_generator/x_electrons_superlattice/sample.yml')
solver = spsc_core.LatticeSymmetrySolver(state)
solver.solve()
state = solver.state
state.export_file('data/solution.yml')


# plt.plot(state.electron_states[0].static_potential)
# plt.plot(state.electron_states[1].static_potential)
# plt.plot(state.static_density * 10 ** (-16))
# plt.show()
# state.export_file('data/RC1402.yml')
# state = spsc_core.StateSimple.import_file('data/solution.yml')
# solver = spsc_core.LatticeXElectronsSymmetrySolver(state)
# solver.solve()
# state.export_file('data/solution.yml')
# state = spsc_data_generator.superlattice_well_sloped(8, spsc_data.LengthValue(26, "nm"),
#                                                      spsc_data.LengthValue(2.3, "nm"), spsc_data.LengthValue(1.1, "nm"),
#                                                      spsc_data.ElectricFieldValue(1.2 * 10 ** -4, "V_per_m"))
# solver = spsc_core.LatticeSlopedSolver(state)
# solver.solve()
# well = {
#     'periods_number': 9,
#     'well_length': 22.0,
#     'lattice_well_length': 2.3,
#     'lattice_barrier_length': 1.4,
#     'delta_layer_period': 4,
#     'delta_layer_density': 2.0 * 10 ** 16,
#     'lattice_amplitude': 1.0,
#     'mass': 0.068,
#     'density': 10.0 * 10 ** 15,
#     'x_lattice_amplitude': 0.295,
#     'x_lattice_offset': 0.46,
#     'mass_x': 1.2
# }
#
# f = open('data/data_generator/x_electrons_superlattice/sample.yml', "w")
# yaml.dump(well, f)
# f.close()
