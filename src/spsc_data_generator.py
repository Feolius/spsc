import spsc_data
import spsc_core
import numpy as np
import spsc_constants as constants


def single_well(length_nm):
    state = spsc_core.StateSimple()
    length = spsc_data.LengthValue(length_nm * 3, "nm")
    granularity = length_nm * 10
    potential = spsc_data.Potential(np.ones((granularity,), "float64"), "eV")
    potential.append(spsc_data.Potential(np.zeros((granularity,), "float64"), "eV"))
    potential.append(spsc_data.Potential(np.ones((granularity + 1,), "float64"), "eV"))
    electron_state = spsc_core.ElectronStateSimple()
    electron_state.static_potential = potential
    electron_state.mass = spsc_data.MassValue(0.068 * constants.m_e)
    state.electron_states = [electron_state]
    state.length = length
    return state


def single_well_sloped(length_nm, slope_V_per_m):
    e = spsc_data.ChargeValue(constants.e)
    length = spsc_data.LengthValue(length_nm * 3, "nm")
    length.convert_to("m")
    final_value_eV = slope_V_per_m * length.value * e.value
    state = single_well(length_nm)
    potential_arr = state.electron_states[0].static_potential.value
    slope_arr = np.zeros((len(potential_arr),))
    dE = final_value_eV / (len(slope_arr) - 1)
    for i in range(len(slope_arr)):
        slope_arr[i] = dE * i
    slope = spsc_data.Potential(slope_arr, "eV")
    potential_arr = potential_arr + slope
    well_min_index = potential_arr.argmin()
    well_min_value = potential_arr[well_min_index]
    potential_arr = potential_arr - np.full((len(potential_arr),), well_min_value)
    state.electron_states[0].static_potential = spsc_data.Potential(potential_arr)
    return state


