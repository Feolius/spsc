import spsc_data
import spsc_core
import numpy as np
import spsc_constants as constants
import yaml
import copy

DOTS_PER_NM = 10


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
    lattice_length = spsc_data.LengthValue(0, "nm")
    for i in range(periods_num):
        lattice_length += lattice_well_length + lattice_barrier_length
    length += lattice_barrier_length * 2 + well_length + lattice_length * 2
    potential = spsc_data.Potential(np.ones((int(lattice_barrier_length.value * DOTS_PER_NM),), "float64"), "eV")
    next_index = lattice_barrier_length.value * DOTS_PER_NM
    for i in range(periods_num):
        potential.append(
            spsc_data.Potential(np.zeros((int(lattice_well_length.value * DOTS_PER_NM),), "float64"), "eV"))
        potential.append(
            spsc_data.Potential(np.ones((int(lattice_barrier_length.value * DOTS_PER_NM),), "float64"), "eV"))
        next_index += lattice_well_length.value * DOTS_PER_NM + lattice_barrier_length.value * DOTS_PER_NM
    meta_info = {
        "well_start": int(next_index),
        "lattice_bar_width": int(lattice_barrier_length.value * DOTS_PER_NM),
        "lattice_well_width": int(lattice_well_length.value * DOTS_PER_NM)
    }
    potential.append(spsc_data.Potential(np.zeros((int(well_length.value * DOTS_PER_NM),), "float64"), "eV"))
    next_index += well_length.value * DOTS_PER_NM
    meta_info["well_end"] = int(next_index)
    empty_dots = int(length.value * DOTS_PER_NM) + 1 - len(potential)
    potential.append(spsc_data.Potential(np.zeros((empty_dots,), "float64"), "eV"))
    potential.mirror()
    potential.meta_info = meta_info
    electron_state.static_potential = potential
    electron_state.mass = spsc_data.MassValue(constants.m_e)
    electron_state.sum_density = spsc_data.DensityValue(0, "m^-2")
    state.electron_states = [electron_state]
    state.length = length
    state.static_density = spsc_data.Density(np.zeros((len(potential, )), "float64"), "m^-2")
    state.density_potential = spsc_data.Potential(np.zeros((len(potential, )), "float64"))
    return state


def superlattice_well_sloped(periods_num, well_length, lattice_well_length, lattice_barrier_length, slope):
    state = superlattice_well(periods_num, well_length, lattice_well_length, lattice_barrier_length)
    slope.convert_to(slope.units_default)
    state.length.convert_to(state.length.units_default)
    energy_diff = slope.value * state.length.value * constants.e
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
    meta_info = state.electron_states[0].static_potential.meta_info
    state.electron_states[0].static_potential = spsc_data.Potential(potential_arr, "eV")
    state.electron_states[0].static_potential.meta_info = meta_info
    return state


def simple_superlattice(file):
    f = open(file, 'r')
    data = yaml.load(f)
    f.close()

    well_length = spsc_data.LengthValue(data['well_length'], 'nm')
    lattice_well_length = spsc_data.LengthValue(data['lattice_well_length'], 'nm')
    lattice_barrier_length = spsc_data.LengthValue(data['lattice_barrier_length'], 'nm')

    state = superlattice_well(data['periods_number'], well_length, lattice_well_length, lattice_barrier_length)
    electron_state = state.electron_states[0]
    electron_state.mass = electron_state.mass * data['mass']
    electron_state.sum_density = spsc_data.DensityValue(data['density'], "m^-2")
    electron_state.static_potential.value = electron_state.static_potential.value * data['lattice_amplitude']

    lattice_well_length.convert_to("nm")
    lattice_barrier_length.convert_to("nm")
    delta_layer_index = (data['delta_layer_period'] - 1) * (
        lattice_well_length.value + lattice_barrier_length.value) * DOTS_PER_NM + \
        lattice_well_length.value * DOTS_PER_NM / 2 + lattice_barrier_length.value * DOTS_PER_NM
    delta_layer_index = int(delta_layer_index)
    state.static_density.convert_to('m^-2')
    state.static_density[delta_layer_index] = data['delta_layer_density']
    state.static_density.mirror()
    state.electron_states[0].static_potential.meta_info['delta_layer_index'] = delta_layer_index
    return state


def x_electrons_superlattice(file):
    f = open(file, 'r')
    data = yaml.load(f)
    f.close()

    well_length = spsc_data.LengthValue(data['well_length'], 'nm')
    lattice_well_length = spsc_data.LengthValue(data['lattice_well_length'], 'nm')
    lattice_barrier_length = spsc_data.LengthValue(data['lattice_barrier_length'], 'nm')

    state = superlattice_well(data['periods_number'], well_length, lattice_well_length, lattice_barrier_length)
    x_electron_state = state.electron_states[0]
    x_electron_state.mass = x_electron_state.mass * data['mass_x']
    x_electron_state.sum_density = spsc_data.DensityValue(0, "m^-2")
    x_electron_state.static_potential.value = x_electron_state.static_potential.value - 1
    x_electron_state.static_potential.value = x_electron_state.static_potential.value * (-data['x_lattice_amplitude'])
    x_electron_state.static_potential.value = x_electron_state.static_potential.value + data['x_lattice_offset']

    state = simple_superlattice(file)

    # Denote start and end indexes for X-electron solutions
    delta_layer_index = state.electron_states[0].static_potential.meta_info['delta_layer_index']
    delta_layer_offset = 2 * (lattice_well_length.value + lattice_barrier_length.value) * DOTS_PER_NM
    x_solution_start_index = delta_layer_index - delta_layer_offset
    if x_solution_start_index < 0:
        x_solution_start_index = 0
    x_solution_end_index = delta_layer_index + delta_layer_offset
    if x_solution_end_index > (len(x_electron_state.static_potential) - 1):
        x_solution_end_index = len(x_electron_state.static_potential) - 1
    x_electron_state.static_potential.meta_info['x_solution_start'] = int(x_solution_start_index)
    x_electron_state.static_potential.meta_info['x_solution_end'] = int(x_solution_end_index)
    x_electron_state.static_potential.meta_info['middle_index'] = int(data['middle_index'])

    state.electron_states.append(x_electron_state)
    return state


def simple_superlattice_diff_mass(file):
    f = open(file, 'r')
    data = yaml.load(f)
    f.close()

    well_length = spsc_data.LengthValue(data['well_length'], 'nm')
    lattice_well_length = spsc_data.LengthValue(data['lattice_well_length'], 'nm')
    lattice_barrier_length = spsc_data.LengthValue(data['lattice_barrier_length'], 'nm')

    state = superlattice_well(data['periods_number'], well_length, lattice_well_length, lattice_barrier_length)
    electron_state = state.electron_states[0]

    mass = spsc_data.MassArray(data['barrier_mass'] * np.ones((int(lattice_barrier_length.value * DOTS_PER_NM),), "float64"))
    for i in range(data['periods_number']):
        mass.append(
            spsc_data.MassArray(data['well_mass'] * np.ones((int(lattice_well_length.value * DOTS_PER_NM),), "float64")))
        mass.append(
            spsc_data.MassArray(data['barrier_mass'] * np.ones((int(lattice_barrier_length.value * DOTS_PER_NM),), "float64")))
    mass.append(spsc_data.MassArray(data['well_mass'] * np.ones((int(well_length.value * DOTS_PER_NM),), "float64")))
    empty_dots = len(electron_state.static_potential) - len(mass)
    mass.append(spsc_data.MassArray(np.zeros((empty_dots,), "float64")))
    mass = mass * constants.m_e
    mass.mirror()
    electron_state.mass = mass

    electron_state.sum_density = spsc_data.DensityValue(data['density'], "m^-2")
    electron_state.static_potential.value = electron_state.static_potential.value * data['lattice_amplitude']

    lattice_well_length.convert_to("nm")
    lattice_barrier_length.convert_to("nm")
    delta_layer_index = (data['delta_layer_period'] - 1) * (
        lattice_well_length.value + lattice_barrier_length.value) * DOTS_PER_NM + \
                        lattice_well_length.value * DOTS_PER_NM / 2 + lattice_barrier_length.value * DOTS_PER_NM
    delta_layer_index = int(delta_layer_index)
    state.static_density.convert_to('m^-2')
    state.static_density[delta_layer_index] = data['delta_layer_density']
    state.static_density.mirror()
    state.electron_states[0].static_potential.meta_info['delta_layer_index'] = delta_layer_index
    return state


def x_electrons_superlattice_diff_mass(file):
    f = open(file, 'r')
    data = yaml.load(f)
    f.close()

    well_length = spsc_data.LengthValue(data['well_length'], 'nm')
    lattice_well_length = spsc_data.LengthValue(data['lattice_well_length'], 'nm')
    lattice_barrier_length = spsc_data.LengthValue(data['lattice_barrier_length'], 'nm')

    state = superlattice_well(data['periods_number'], well_length, lattice_well_length, lattice_barrier_length)
    x_electron_state = state.electron_states[0]
    x_electron_state.mass = x_electron_state.mass * data['well_mass_x']
    x_electron_state.sum_density = spsc_data.DensityValue(0, "m^-2")
    x_electron_state.static_potential.value = x_electron_state.static_potential.value - 1
    x_electron_state.static_potential.value = x_electron_state.static_potential.value * (-data['x_lattice_amplitude'])
    x_electron_state.static_potential.value = x_electron_state.static_potential.value + data['x_lattice_offset']

    state = simple_superlattice_diff_mass(file)

    # Denote start and end indexes for X-electron solutions
    delta_layer_index = state.electron_states[0].static_potential.meta_info['delta_layer_index']
    delta_layer_offset = 2 * (lattice_well_length.value + lattice_barrier_length.value) * DOTS_PER_NM
    x_solution_start_index = delta_layer_index - delta_layer_offset
    if x_solution_start_index < 0:
        x_solution_start_index = 0
    x_solution_end_index = delta_layer_index + delta_layer_offset
    if x_solution_end_index > (len(x_electron_state.static_potential) - 1):
        x_solution_end_index = len(x_electron_state.static_potential) - 1
    x_electron_state.static_potential.meta_info['x_solution_start'] = int(x_solution_start_index)
    x_electron_state.static_potential.meta_info['x_solution_end'] = int(x_solution_end_index)
    x_electron_state.static_potential.meta_info['middle_index'] = int(data['middle_index'])

    state.electron_states.append(x_electron_state)
    return state
