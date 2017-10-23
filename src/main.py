import spsc_data_generator
import spsc_core
import spsc_data
import spsc_constants
import numpy as np
import matplotlib.pyplot as plt

# state = spsc_data_generator.simple_superlattice_diff_mass('data/data_generator/x_electrons_superlattice_diff_mass/sample.yml')
# electron_state = state.electron_states[0]
# plt.plot(electron_state.mass)
# electron_state.static_potential.value = electron_state.static_potential.value * electron_state.mass.value[0]
# plt.plot(electron_state.static_potential)
# plt.show()

# state = spsc_data_generator.x_electrons_superlattice_diff_mass('data/data_generator/x_electrons_superlattice_diff_mass/sample.yml')
# solver = spsc_core.LatticeXElectronsSymmetrySolver(state)
# solver.solve()
# state = solver.state
# state.export_file('data/solutions/x_electrons_superlattice/RC1075.yml')
#
state = spsc_core.StateSimple.import_file('data/solutions/x_electrons_superlattice/RC1075.yml')
electron_state_1 = state.electron_states[0]
potential_1 = electron_state_1.static_potential + state.density_potential
potential_1.convert_to("eV")
E1 = electron_state_1.energy_levels[0]
E1.convert_to('eV')
E2 = electron_state_1.energy_levels[1]
E2.convert_to('eV')
plt.plot(potential_1)
wf_1 = electron_state_1.wave_functions[0]
plt.plot(wf_1)
wf_2 = electron_state_1.wave_functions[1]
plt.plot(wf_2)

electron_state_2 = state.electron_states[1]
potential_2 = electron_state_2.static_potential + state.density_potential
potential_2.convert_to('eV')
plt.plot(potential_2)
wf_x = electron_state_2.wave_functions[0]
wf_x.value = wf_x.value * 0.2
Ex = electron_state_2.energy_levels[0]
Ex.convert_to('eV')
plt.plot(wf_x)

electron_state_1.static_potential.convert_to('eV')
plt.plot(E1.value * np.ones((len(potential_1), ), "float64"))
plt.plot(E2.value * np.ones((len(potential_1), ), "float64"))
plt.plot(Ex.value * np.ones((len(potential_1), ), "float64"))

plt.show()

n = spsc_data.DensityValue(6.81 * 10 ** 15, 'm^-2')
n.convert_to(n.units_default)
electron_state_1.mass.convert_to(electron_state_1.mass.units_default)
g = 0.068 * spsc_constants.m_e / (np.pi * spsc_constants.h_plank ** 2)
E1.convert_to(E1.units_default)
E2.convert_to(E2.units_default)
Ef = spsc_data.EnergyValue(0.5*(n.value / g + E1.value + E2.value))
Ef1 = Ef - E1
Ef2 = Ef - E2
Ef.convert_to('eV')
E1.convert_to('eV')
E2.convert_to('eV')
Ef1.convert_to('eV')
Ef2.convert_to('eV')
print "E1: ", E1
print "E2: ", E2
print "dE: ", (E2 - E1)
print "Ef1: ", Ef1
print "Ef2: ", Ef2
print "Ex: ", Ex

# E1:  0.209943854689 eV
# E2:  0.220043854689 eV
# dE:  0.0101 eV
# Ef1:  0.0170370147962 eV
# Ef2:  0.00693701479618 eV
# Ex:  0.1915 eV

# import start_scripts.RC1075