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
electron_state_1 = state.electron_states[0]
potential_1 = electron_state_1.static_potential + state.density_potential
potential_1.convert_to("eV")
electron_state_2 = state.electron_states[1]
potential_2 = electron_state_2.static_potential + state.density_potential
potential_2.convert_to("eV")
length = state.length
length.convert_to("nm")
x = np.linspace(0, length.value, len(potential_1))
ax = plt.subplot(111)
ax.plot(x, potential_1, label=r'$\Gamma$')
ax.plot(x, potential_2, label="X")
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Only show ticks on the left and bottom spines
ax.yaxis.set_ticks_position('left')
ax.xaxis.set_ticks_position('bottom')
ax.set_xlabel('Z, nm', fontsize=16)
ax.set_ylabel('V, eV', fontsize=16)
ax.legend((r'$\Gamma$', 'X'), fontsize=11)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(16)
# plt.savefig("GaAs_AlAs_init_potential.png", dpi=400, bbox_inches='tight')
plt.show()
solver = spsc_solver.LatticeXElectronsSymmetrySolver(state)
solver.solve()
state = solver.state
# state.export_file('data/solutions/x_electrons_superlattice_diff_mass/disser.yml')

# state = spsc_core.StateSimple.import_file('data/solutions/x_electrons_superlattice_diff_mass/disser.yml')
length = state.length
length.convert_to("nm")
electron_state_1 = state.electron_states[0]
potential_1 = electron_state_1.static_potential + state.density_potential
potential_1.convert_to("eV")
E1 = electron_state_1.energy_levels[0]
E1.convert_to('eV')
wf_1 = electron_state_1.wave_functions[0]
wf_1.value = wf_1.value * 0.3

electron_state_2 = state.electron_states[1]
potential_2 = electron_state_2.static_potential + state.density_potential
potential_2.convert_to('eV')
wf_x = electron_state_2.wave_functions[0]
wf_x.value = wf_x.value * 0.1
Ex = electron_state_2.energy_levels[0]
Ex.convert_to('eV')

electron_state_1.static_potential.convert_to('eV')

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

x = np.linspace(0, length.value, len(potential_1))
ax = plt.subplot(111)
# Hide the right and top spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.set_xlabel('Z, nm', fontsize=16)
ax.set_ylabel('V, eV', fontsize=16)
ax.plot(x, potential_1, label=r'$\Gamma$')
ax.plot(x, potential_2, label='X')
ax.plot(x, wf_1, label=r'$\psi_\Gamma$')
ax.plot(x, wf_x, label=r'$\psi_X$')
ax.plot(x, E1.value * np.ones((len(potential_1), ), "float64"), label=r'$E_\Gamma$')
ax.plot(x, Ex.value * np.ones((len(potential_1), ), "float64"), label=r'$E_x$')
ax.plot(x, Ef.value * np.ones((len(potential_1), ), "float64"), label=r'$E_f$', ls='-.', c='black')
ax.legend(fontsize=14)
for tick in ax.xaxis.get_major_ticks():
    tick.label.set_fontsize(16)
for tick in ax.yaxis.get_major_ticks():
    tick.label.set_fontsize(16)
plt.savefig("GaAs_AlAs_solution.png", dpi=400, bbox_inches='tight')
plt.show()

print "E1: ", E1
print "Ef1: ", Ef1
print "Ex: ", Ex

