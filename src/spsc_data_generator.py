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
    e = slope_arr[9] - slope_arr[0]
    potential_arr = potential_arr + slope_arr
    well_min_index = potential_arr.argmin()
    well_min_value = potential_arr[well_min_index]
    e = potential_arr[9] - potential_arr[0]
    potential_arr = potential_arr - np.full((len(potential_arr),), well_min_value)
    e = potential_arr[9] - potential_arr[0]
    state.electron_states[0].static_potential = spsc_data.Potential(potential_arr, "eV")
    return state


