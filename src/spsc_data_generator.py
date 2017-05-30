import spsc_data
import spsc_core
import numpy as np
import spsc_constants as constants
import copy


def single_well(well_length):
    state = spsc_core.StateSimple()
    granularity = well_length.value * 10
    length = well_length * 3
    length.convert_to("nm")
    potential = spsc_data.Potential(np.ones((granularity,), "float64"), "eV")
    potential.append(spsc_data.Potential(np.zeros((granularity,), "float64"), "eV"))
    potential.append(spsc_data.Potential(np.ones((granularity + 1,), "float64"), "eV"))
    electron_state = spsc_core.ElectronStateSimple()
    electron_state.static_potential = potential
    electron_state.mass = spsc_data.MassValue(0.068 * constants.m_e)
    state.electron_states = [electron_state]
    state.length = length
    return state


def single_well_sloped(well_length, slope):
    state = single_well(well_length)
    length = well_length * 3
    slope.convert_to(slope.units_default)
    length.convert_to(length.units_default)
    energy_diff = slope.value * length.value * constants.e
    energy_diff = spsc_data.EnergyValue(energy_diff)
    energy_diff.convert_to("eV")

    state.electron_states[0].static_potential.convert_to("eV")
    potential_arr = state.electron_states[0].static_potential.value
    slope_arr = np.zeros((len(potential_arr),))
    dE = energy_diff.value / (len(slope_arr) - 1)
    for i in range(len(slope_arr)):
        slope_arr[i] = dE * i
    potential_arr = potential_arr + slope_arr
    well_min_index = potential_arr.argmin()
    well_min_value = potential_arr[well_min_index]
    potential_arr = potential_arr - np.full((len(potential_arr),), well_min_value)
    state.electron_states[0].static_potential = spsc_data.Potential(potential_arr, "eV")
    return state


def superlattice_well(periods_num, well_length, lattice_well_length, lattice_barrier_length):
    well_length.convert_to("nm")
    lattice_well_length.convert_to("nm")
    lattice_barrier_length.convert_to("nm")
    state = spsc_core.StateSimple()
    electron_state = spsc_core.ElectronStateSimple()
    length = spsc_data.LengthValue(0, "nm")
    dots_per_nm = 10
    # Let's make a surface layer with well length
    lattice_length = spsc_data.LengthValue(0, "nm")
    for i in range(periods_num):
        lattice_length += lattice_well_length + lattice_barrier_length
    length += well_length * 3 + lattice_length * 2
    potential = spsc_data.Potential(5 * np.ones((well_length.value * dots_per_nm,), "float64"), "eV")
    for i in range(periods_num):
        potential.append(spsc_data.Potential(np.zeros((int(lattice_well_length.value * dots_per_nm),), "float64"), "eV"))
        potential.append(spsc_data.Potential(np.ones((int(lattice_barrier_length.value * dots_per_nm),), "float64"), "eV"))
    potential.append(spsc_data.Potential(np.zeros((int(well_length.value * dots_per_nm),), "float64"), "eV"))
    empty_dots = int(length.value * dots_per_nm) + 1 - len(potential)
    potential.append(spsc_data.Potential(np.zeros((empty_dots,), "float64"), "eV"))
    potential.mirror()
    electron_state.static_potential = potential
    electron_state.mass = spsc_data.MassValue(0.068 * constants.m_e)
    state.electron_states = [electron_state]
    state.length = length
    state.static_density = spsc_data.Density(np.zeros((len(potential,)), "float64"), "m^-2")
    density_index = (well_length.value + lattice_length.value / 2) * dots_per_nm
    state.static_density[int(density_index)] = 1.6 * 10 ** 16
    state.static_density.mirror()
    return state


