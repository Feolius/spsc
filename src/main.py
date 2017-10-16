import spsc_data_generator
import matplotlib.pyplot as plt

state = spsc_data_generator.simple_superlattice_diff_mass('data/data_generator/simple_superlattice_diff_mass/sample.yml')
electron_state = state.electron_states[0]
plt.plot(electron_state.mass)
plt.show()

