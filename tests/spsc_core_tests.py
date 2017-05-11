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


def random_electron_state_simple():
    energy_levels_number = random.randint(1, 5)
    energy_levels = []
    rnd_energy_units = random_units("EnergyValue")
    for i in range(energy_levels_number):
        energy_levels.append(random_value("EnergyValue", rnd_energy_units))
    wave_functions = []
    rnd_granularity = random.randint(10, 20)
    for i in range(energy_levels_number):
        wave_functions.append(random_value_array("WaveFunction", rnd_granularity))
    static_potential = random_value_array("Potential", rnd_granularity, rnd_energy_units)
    density_potential = random_value_array("Potential", rnd_granularity, rnd_energy_units)
    mass = random_value("MassValue")
    electron_state = spsc_core.ElectronStateSimple()
    electron_state.energy_levels = energy_levels
    electron_state.wave_functions = wave_functions
    electron_state.static_potential = static_potential
    electron_state.density_potential = density_potential
    electron_state.mass = mass
    return electron_state


class ElectronStateSimpleOperatorTestCase(unittest.TestCase):

    def test_equals(self):
        electron_state1 = random_electron_state_simple()
        electron_state2 = copy.deepcopy(electron_state1)
        self.assertEqual(electron_state1, electron_state2)

        electron_state2.mass += spsc_data.MassValue(1)
        self.assertNotEqual(electron_state1, electron_state2)
        electron_state2 = copy.deepcopy(electron_state1)

        electron_state2.density_potential += random_value_array("Potential", len(electron_state2.density_potential),
                                                                electron_state2.density_potential.units)
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
        dct_density_potential = spsc_data.Potential.from_dict(dct["density_potential"])
        self.assertEqual(dct_density_potential, electron_state.density_potential)
        dct_mass = spsc_data.MassValue.from_dict(dct["mass"])
        self.assertEqual(dct_mass, electron_state.mass)

    def test_from_dict(self):
        energy_levels_number = random.randint(1, 5)
        rnd_energy_units = random_units("EnergyValue")
        energy_levels = []
        for i in range(energy_levels_number):
            energy_levels.append(random_value("EnergyValue", rnd_energy_units).to_dict())
        rnd_granularity = random.randint(10, 20)
        wave_functions = []
        for i in range(energy_levels_number):
            wave_functions.append(random_value_array("WaveFunction", rnd_granularity).to_dict())
        static_potential = random_value_array("Potential", rnd_granularity, rnd_energy_units).to_dict()
        density_potential = random_value_array("Potential", rnd_granularity, rnd_energy_units).to_dict()
        mass = random_value("MassValue").to_dict()
        electron_state_dict = {
            "energy_levels": energy_levels,
            "wave_functions": wave_functions,
            "static_potential": static_potential,
            "density_potential": density_potential,
            "mass": mass
        }
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
        dict_density_potential = spsc_data.Potential.from_dict(electron_state_dict["density_potential"])
        self.assertEqual(electron_state.density_potential, dict_density_potential)
        dict_static_potential = spsc_data.Potential.from_dict(electron_state_dict["static_potential"])
        self.assertEqual(electron_state.static_potential, dict_static_potential)
        dict_mass = spsc_data.MassValue.from_dict(electron_state_dict["mass"])
        self.assertEqual(electron_state.mass, dict_mass)

    def test_import_export(self):
        electron_state = random_electron_state_simple()
        test_file = "test.yml"
        electron_state.export_file(test_file)
        electron_state_imported = spsc_core.ElectronStateSimple.import_file(test_file)
        self.assertEqual(electron_state, electron_state_imported)
        os.remove(test_file)