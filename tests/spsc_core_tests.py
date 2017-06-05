import unittest
import src.spsc_core as spsc_core
import src.spsc_data as spsc_data
import random
import copy
import os


def random_list(length):
    rnd_lst = []
    for i in range(length):
        rnd_lst.append(random.uniform(0, 100))
    return rnd_lst


def random_units(class_name):
    class_ = getattr(spsc_data, class_name)
    units = class_.get_all_units()
    rnd_unit = units[random.randrange(0, len(units))]
    return rnd_unit


def random_value(class_name, units=None):
    class_ = getattr(spsc_data, class_name)
    if units is None:
        units = random_units(class_name)
    value = class_(random.uniform(1, 100), units)
    return value


def random_value_array(class_name, length, units=None):
    class_ = getattr(spsc_data, class_name)
    if units is None:
        units = random_units(class_name)
    arr = class_(random_list(length), units)
    return arr


def random_electron_state_simple(granularity=None):
    energy_levels_number = random.randint(1, 5)
    energy_levels = []
    rnd_energy_units = random_units("EnergyValue")
    for i in range(energy_levels_number):
        energy_levels.append(random_value("EnergyValue", rnd_energy_units))
    wave_functions = []
    if granularity is None:
        granularity = random.randint(10, 20)
    for i in range(energy_levels_number):
        wave_functions.append(random_value_array("WaveFunction", granularity))
    static_potential = random_value_array("Potential", granularity, rnd_energy_units)
    mass = random_value("MassValue")
    sum_density = random_value("DensityValue")
    electron_state = spsc_core.ElectronStateSimple()
    electron_state.energy_levels = energy_levels
    electron_state.wave_functions = wave_functions
    electron_state.static_potential = static_potential
    electron_state.mass = mass
    electron_state.sum_density = sum_density
    return electron_state


def random_electron_state_simple_dict(granularity=None):
    energy_levels_number = random.randint(1, 5)
    rnd_energy_units = random_units("EnergyValue")
    energy_levels = []
    for i in range(energy_levels_number):
        energy_levels.append(random_value("EnergyValue", rnd_energy_units).to_dict())
    if granularity is None:
        granularity = random.randint(10, 20)
    wave_functions = []
    for i in range(energy_levels_number):
        wave_functions.append(random_value_array("WaveFunction", granularity).to_dict())
    static_potential = random_value_array("Potential", granularity, rnd_energy_units).to_dict()
    mass = random_value("MassValue").to_dict()
    sum_density = random_value("DensityValue").to_dict()
    dct = {
        "energy_levels": energy_levels,
        "wave_functions": wave_functions,
        "static_potential": static_potential,
        "mass": mass,
        "sum_density": sum_density
    }
    return dct


def random_state_simple():
    state = spsc_core.StateSimple()
    rnd_electron_states_num = random.randrange(1, 5)
    rnd_granularity = random.randint(10, 20)
    electron_states = []
    for i in range(rnd_electron_states_num):
        electron_states.append(random_electron_state_simple(rnd_granularity))
    state.electron_states = electron_states
    state.static_density = random_value_array("Density", rnd_granularity, random_units("Density"))
    state.density_potential = random_value_array("Potential", rnd_granularity, random_units("EnergyValue"))
    state.length = random_value("LengthValue", random_units("LengthValue"))
    return state


def random_state_simple_dict():
    electron_states = []
    rnd_electron_states_num = random.randrange(1, 5)
    rnd_granularity = random.randint(10, 20)
    for i in range(rnd_electron_states_num):
        electron_states.append(random_electron_state_simple(rnd_granularity).to_dict())
    dct = {
        "electron_states": electron_states,
        "static_density": random_value_array("Density", rnd_granularity, random_units("Density")).to_dict(),
        "density_potential": random_value_array("Potential", rnd_granularity, random_units("EnergyValue")).to_dict(),
        "length": random_value("LengthValue", random_units("LengthValue")).to_dict()
    }
    return dct


class ElectronStateSimpleHelpersTestCase(unittest.TestCase):

    def test_is_granularity_consistent(self):
        electron_state_original = random_electron_state_simple()
        self.assertTrue(electron_state_original._is_granularity_consistent())
        rnd_int = random.randrange(0, len(electron_state_original.static_potential))

        electron_state = copy.deepcopy(electron_state_original)
        del electron_state.static_potential[rnd_int]
        self.assertFalse(electron_state._is_granularity_consistent())

        for i in range(len(electron_state_original.wave_functions)):
            electron_state = copy.deepcopy(electron_state_original)
            del electron_state.wave_functions[i][rnd_int]
            self.assertFalse(electron_state._is_granularity_consistent())

    def test_is_wave_functions_consistent(self):
        electron_state_original = random_electron_state_simple()
        self.assertTrue(electron_state_original._is_wave_functions_consistent())
        rnd_int = random.randrange(0, len(electron_state_original.wave_functions))

        electron_state = copy.deepcopy(electron_state_original)
        del electron_state.wave_functions[rnd_int]
        self.assertFalse(electron_state._is_wave_functions_consistent())

        del electron_state_original.energy_levels[rnd_int]
        self.assertFalse(electron_state_original._is_wave_functions_consistent())

    def test_get_granularity(self):
        electron_state = spsc_core.ElectronStateSimple()
        self.assertEqual(electron_state.get_granularity(), 0)
        electron_state = random_electron_state_simple()
        self.assertEqual(electron_state.get_granularity(), len(electron_state.static_potential))
        rnd_int = random.randrange(0, len(electron_state.static_potential))
        del electron_state.static_potential[rnd_int]
        with self.assertRaises(StandardError) as context:
            electron_state.get_granularity()


class ElectronStateSimpleOperatorTestCase(unittest.TestCase):

    def test_equals(self):
        electron_state1 = random_electron_state_simple()
        electron_state2 = copy.deepcopy(electron_state1)
        self.assertEqual(electron_state1, electron_state2)

        electron_state2.mass += spsc_data.MassValue(1)
        self.assertNotEqual(electron_state1, electron_state2)
        electron_state2 = copy.deepcopy(electron_state1)

        electron_state2.static_potential += random_value_array("Potential", len(electron_state2.static_potential),
                                                               electron_state2.static_potential.units)
        self.assertNotEqual(electron_state1, electron_state2)
        electron_state2 = copy.deepcopy(electron_state1)

        rnd_int = random.randrange(0, len(electron_state1.energy_levels))

        electron_state2.energy_levels[rnd_int] += random_value("EnergyValue", electron_state2.energy_levels[rnd_int].units)
        self.assertNotEqual(electron_state1, electron_state2)
        del electron_state2.energy_levels[rnd_int]
        self.assertNotEqual(electron_state1, electron_state2)
        electron_state2 = copy.deepcopy(electron_state1)

        electron_state2.wave_functions[rnd_int] += random_value_array("WaveFunction", len(electron_state2.wave_functions[rnd_int]))
        self.assertNotEqual(electron_state1, electron_state2)


class ElectronStateSimpleIOTestCase(unittest.TestCase):

    def test_to_dict(self):
        electron_state = random_electron_state_simple()
        dct = electron_state.to_dict()
        self.assertEqual(len(dct["energy_levels"]), len(electron_state.energy_levels))
        for i in range(len(dct["energy_levels"])):
            dct_energy_level = spsc_data.EnergyValue.from_dict(dct["energy_levels"][i])
            self.assertEqual(dct_energy_level, electron_state.energy_levels[i])
        self.assertEqual(len(dct["wave_functions"]), len(electron_state.wave_functions))
        for i in range(len(dct["wave_functions"])):
            dct_wave_function = spsc_data.WaveFunction.from_dict(dct["wave_functions"][i])
            self.assertEqual(dct_wave_function, electron_state.wave_functions[i])
        self.assertEqual(len(dct["energy_levels"]), len(dct["wave_functions"]))
        dct_static_potential = spsc_data.Potential.from_dict(dct["static_potential"])
        self.assertEqual(dct_static_potential, electron_state.static_potential)
        dct_mass = spsc_data.MassValue.from_dict(dct["mass"])
        self.assertEqual(dct_mass, electron_state.mass)
        dct_sum_density = spsc_data.DensityValue.from_dict(dct["sum_density"])
        self.assertEqual(dct_sum_density, electron_state.sum_density)

    def test_from_dict(self):
        electron_state_dict = random_electron_state_simple_dict()
        electron_state = spsc_core.ElectronStateSimple.from_dict(electron_state_dict)
        self.assertEqual(len(electron_state.wave_functions), len(electron_state_dict["wave_functions"]))
        for i in range(len(electron_state.wave_functions)):
            dict_wave_function = spsc_data.WaveFunction.from_dict(electron_state_dict["wave_functions"][i])
            self.assertEqual(electron_state.wave_functions[i], dict_wave_function)
        self.assertEqual(len(electron_state.energy_levels), len(electron_state_dict["energy_levels"]))
        for i in range(len(electron_state.energy_levels)):
            dict_energy_level = spsc_data.EnergyValue.from_dict(electron_state_dict["energy_levels"][i])
            self.assertEqual(electron_state.energy_levels[i], dict_energy_level)
        self.assertEqual(len(electron_state.wave_functions), len(electron_state.energy_levels))
        dict_static_potential = spsc_data.Potential.from_dict(electron_state_dict["static_potential"])
        self.assertEqual(electron_state.static_potential, dict_static_potential)
        dict_mass = spsc_data.MassValue.from_dict(electron_state_dict["mass"])
        self.assertEqual(electron_state.mass, dict_mass)
        dict_sum_density = spsc_data.DensityValue.from_dict(electron_state_dict["sum_density"])
        self.assertEqual(electron_state.sum_density, dict_sum_density)

    def test_import_export(self):
        electron_state = random_electron_state_simple()
        test_file = "test.yml"
        electron_state.export_file(test_file)
        electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)
        self.assertEqual(electron_state, electron_state_imported)
        os.remove(test_file)

    def test_import_export_incorrect(self):
        electron_state = random_electron_state_simple()
        rnd_int = random.randrange(0, len(electron_state.energy_levels))
        del electron_state.energy_levels[rnd_int]
        test_file = "test.yml"
        electron_state.export_file(test_file)
        with self.assertRaises(ImportError) as context:
            electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)

        electron_state = random_electron_state_simple()
        rnd_int = random.randrange(0, len(electron_state.energy_levels))
        del electron_state.wave_functions[rnd_int]
        electron_state.export_file(test_file)
        with self.assertRaises(ImportError) as context:
            electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)

        electron_state_original = random_electron_state_simple()
        rnd_int = random.randrange(0, electron_state_original.get_granularity())

        electron_state = copy.deepcopy(electron_state_original)
        del electron_state.static_potential[rnd_int]
        electron_state.export_file(test_file)
        with self.assertRaises(ImportError) as context:
            electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)

        electron_state = copy.deepcopy(electron_state_original)
        rnd_wave_func_index = random.randrange(0, len(electron_state.wave_functions))
        del electron_state.wave_functions[rnd_wave_func_index][rnd_int]
        electron_state.export_file(test_file)
        with self.assertRaises(ImportError) as context:
            electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)

        electron_state = copy.deepcopy(electron_state_original)
        del electron_state.energy_levels[rnd_wave_func_index]
        electron_state.export_file(test_file)
        with self.assertRaises(ImportError) as context:
            electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)


class StateSimpleHelpersTestCase(unittest.TestCase):

    def test_is_granularity_consistent(self):
        state = random_state_simple()
        self.assertTrue(state._is_granularity_consistent())
        granularity = state.electron_states[0].get_granularity()
        rnd_electron_state = random_electron_state_simple(granularity + 1)
        state.electron_states.append(rnd_electron_state)
        self.assertFalse(state._is_granularity_consistent())

    def test_get_granularity(self):
        state = spsc_core.StateSimple()
        self.assertEqual(state.get_granularity(), 0)
        state = random_state_simple()
        granularity = state.electron_states[0].get_granularity()
        self.assertEqual(state.get_granularity(), granularity)
        rnd_electron_state = random_electron_state_simple(granularity + 1)
        state.electron_states.append(rnd_electron_state)
        with self.assertRaises(StandardError) as context:
            state.get_granularity()


class StateSimpleOperatorTestCase(unittest.TestCase):

    def test_equals(self):
        state1 = random_state_simple()
        state2 = copy.deepcopy(state1)
        self.assertEqual(state1, state2)

        granularity = state2.electron_states[0].get_granularity()
        rnd_electron_state = random_electron_state_simple(granularity + 1)
        rnd_index = random.randrange(0, len(state2.electron_states))
        state2.electron_states[rnd_index] = rnd_electron_state
        self.assertNotEqual(state1, state2)

        state2 = copy.deepcopy(state1)
        state2.density_potential += random_value_array("Potential", len(state2.density_potential), state2.density_potential.units)
        self.assertNotEqual(state1, state2)

        state2 = copy.deepcopy(state1)
        state2.static_density += random_value_array("Density", len(state2.static_density))
        self.assertNotEqual(state1, state2)

        state2 = copy.deepcopy(state1)
        state2.length += random_value("LengthValue")
        self.assertNotEqual(state1, state2)


class StateSimpleIOTestCase(unittest.TestCase):

    def test_to_dict(self):
        state = random_state_simple()
        state_dct = state.to_dict()
        self.assertEqual(len(state_dct["electron_states"]), len(state.electron_states))
        for i in range(len(state.electron_states)):
            electron_state_dct = spsc_core.ElectronStateSimple.from_dict(state_dct["electron_states"][i])
            self.assertEqual(electron_state_dct, state.electron_states[i])
        density_potential_dct = spsc_data.Potential.from_dict(state_dct["density_potential"])
        self.assertEqual(state.density_potential, density_potential_dct)
        static_density_dct = spsc_data.Density.from_dict(state_dct["static_density"])
        self.assertEqual(state.static_density, static_density_dct)
        length_dct = spsc_data.LengthValue.from_dict(state_dct["length"])
        self.assertEqual(length_dct, state.length)

    def test_from_dict(self):
        state_dict = random_state_simple_dict()
        state = spsc_core.StateSimple.from_dict(state_dict)
        self.assertEqual(len(state.electron_states), len(state_dict["electron_states"]))
        for i in range(len(state_dict["electron_states"])):
            dict_electron_state = spsc_core.ElectronStateSimple.from_dict(state_dict["electron_states"][i])
            self.assertEqual(dict_electron_state, state.electron_states[i])
        dict_density_potential = spsc_data.Potential.from_dict(state_dict["density_potential"])
        self.assertEqual(dict_density_potential, state.density_potential)
        dict_static_density = spsc_data.Density.from_dict(state_dict["static_density"])
        self.assertEqual(dict_static_density, state.static_density)
        dict_length = spsc_data.LengthValue.from_dict(state_dict["length"])
        self.assertEqual(dict_length, state.length)

    def test_from_dict_incorrect(self):
        state_dict = random_state_simple_dict()
        electron_state_first = spsc_core.ElectronStateSimple.from_dict(state_dict["electron_states"][0])
        granularity = electron_state_first.get_granularity()
        different_granularity = granularity + random.randint(1, 10)
        new_electron_state = random_electron_state_simple_dict(different_granularity)
        state_dict["electron_states"].append(new_electron_state)
        with self.assertRaises(ImportError) as context:
            state = spsc_core.StateSimple.from_dict(state_dict)

    def test_import_export(self):
        state = random_state_simple()
        test_file = "test.yml"
        state.export_file(test_file)
        state_imported = spsc_core.StateSimple.import_file(test_file)
        self.assertEqual(state, state_imported)

