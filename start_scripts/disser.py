import scipy.special as sp
import scipy.linalg as la
import numpy as np
import matplotlib.pyplot as plt
import src.spsc_data_generator as spsc_data_generator
import src.spsc_data as spsc_data
import src.spsc_constants as spsc_constants
import src.spsc_core as spsc_core
import src.spsc_solver as spsc_solver
import yaml

state = spsc_data_generator.x_electrons_superlattice_diff_mass('data/data_generator/x_electrons_superlattice_diff_mass/disser.yml')
# electron_state_1 = state.electron_states[0]
# potential_1 = electron_state_1.static_potential + state.density_potential
# potential_1.convert_to("eV")
# electron_state_2 = state.electron_states[1]
# potential_2 = electron_state_2.static_potential + state.density_potential
# potential_2.convert_to("eV")
# plt.plot(potential_1)
# plt.plot(potential_2)
# plt.show()
solver = spsc_solver.LatticeXElectronsSymmetrySolver(state)
solver.solve()
state = solver.state
# state.export_file('data/solutions/x_electrons_superlattice_diff_mass/disser.yml')

# state = spsc_core.StateSimple.import_file('data/solutions/x_electrons_superlattice_diff_mass/disser.yml')
electron_state_1 = state.electron_states[0]
potential_1 = electron_state_1.static_potential + state.density_potential
potential_1.convert_to("eV")
E1 = electron_state_1.energy_levels[0]
E1.convert_to('eV')
plt.plot(potential_1)
wf_1 = electron_state_1.wave_functions[0]
wf_1.value = wf_1.value * 0.3
plt.plot(wf_1)

electron_state_2 = state.electron_states[1]
potential_2 = electron_state_2.static_potential + state.density_potential
potential_2.convert_to('eV')
plt.plot(potential_2)
wf_x = electron_state_2.wave_functions[0]
wf_x.value = wf_x.value * 0.1
Ex = electron_state_2.energy_levels[0]
Ex.convert_to('eV')
plt.plot(wf_x)

electron_state_1.static_potential.convert_to('eV')
plt.plot(E1.value * np.ones((len(potential_1), ), "float64"))
plt.plot(Ex.value * np.ones((len(potential_1), ), "float64"))

plt.show()

n = spsc_data.DensityValue(6.81 * 10 ** 15, 'm^-2')
n.convert_to(n.units_default)
electron_state_1.mass.convert_to(electron_state_1.mass.units_default)
well_start = electron_state_1.static_potential.meta_info['well_start']
g = electron_state_1.mass.value[well_start] / (np.pi * spsc_constants.h_plank ** 2)
E1.convert_to(E1.units_default)
Ef = spsc_data.EnergyValue(n.value / g + E1.value)
Ef1 = Ef - E1
Ef.convert_to('eV')
E1.convert_to('eV')
Ef1.convert_to('eV')
print "E1: ", E1
print "Ef1: ", Ef1
print "Ex: ", Ex

