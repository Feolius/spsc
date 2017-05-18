import spsc_data
import spsc_core
import numpy as np
import scipy.constants as constants


def single_well():
    state = spsc_core.StateSimple()
    length = spsc_data.LengthValue(26 * 3, "nm")
    potential = spsc_data.Potential(np.ones((260,), "float64"), "eV")
    potential.append(spsc_data.Potential(np.zeros((260,), "float64"), "eV"))
    potential.append(spsc_data.Potential(np.ones((261,), "float64"), "eV"))
    electron_state = spsc_core.ElectronStateSimple()
    electron_state.static_potential = potential
    electron_state.mass = spsc_data.MassValue(0.068 * constants.m_e, "kg")
    state.electron_states = [electron_state]
    state.length = length
    return state
