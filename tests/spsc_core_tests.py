import unittest
import src.spsc_core as spsc_core
import src.spsc_data as spsc_data
import random
import os


class ElectronStateSimpleIOTestCase(unittest.TestCase):

    @staticmethod
    def random_list(length):
        rnd_lst = []
        for i in range(length):
            rnd_lst.append(random.uniform(0, 100))
        return rnd_lst

    @staticmethod
    def random_units(class_name):
        class_ = getattr(spsc_data, class_name)
        units = class_.get_all_units()
        rnd_unit = units[random.randrange(0, len(units))]
        return rnd_unit

    @staticmethod
    def random_value(class_name, units=None):
        class_ = getattr(spsc_data, class_name)
        if units is None:
            units = ElectronStateSimpleIOTestCase.random_units(class_name)
        value = class_(random.uniform(0, 100), units)
        return value

    @staticmethod
    def random_value_array(class_name, length, units=None):
        class_ = getattr(spsc_data, class_name)
        if units is None:
            units = ElectronStateSimpleIOTestCase.random_units(class_name)
        arr = class_(ElectronStateSimpleIOTestCase.random_list(length), units)
        return arr

    def test_to_dict(self):
        energy_levels_number = random.randint(1, 5)
        energy_levels = []
        rnd_energy_units = ElectronStateSimpleIOTestCase.random_units("EnergyValue")
        for i in range(energy_levels_number):
            energy_levels.append(ElectronStateSimpleIOTestCase.random_value("EnergyValue", rnd_energy_units))
        wave_functions = []
        rnd_granularity = random.randint(10, 20)
        for i in range(energy_levels_number):
            wave_functions.append(ElectronStateSimpleIOTestCase.random_value_array("WaveFunction", rnd_granularity))
        static_potential = ElectronStateSimpleIOTestCase.random_value_array("Potential", rnd_granularity,
                                                                            rnd_energy_units)
        density_potential = ElectronStateSimpleIOTestCase.random_value_array("Potential", rnd_granularity,
                                                                             rnd_energy_units)
        mass = ElectronStateSimpleIOTestCase.random_value("MassValue")
        electron_state = spsc_core.ElectronStateSimple()
        electron_state.energy_levels = energy_levels
        electron_state.wave_functions = wave_functions
        electron_state.static_potential = static_potential
        electron_state.density_potential = density_potential
        electron_state.mass = mass
        dct = electron_state.to_dict()
        self.assertEqual(len(dct["energy_levels"]), len(energy_levels))
        for i in range(len(dct["energy_levels"])):
            dct_energy_level = spsc_data.EnergyValue.from_dict(dct["energy_levels"][i])
            self.assertEqual(dct_energy_level, energy_levels[i])
        self.assertEqual(len(dct["wave_functions"]), len(wave_functions))
        for i in range(len(dct["wave_functions"])):
            dct_wave_function = spsc_data.WaveFunction.from_dict(dct["wave_functions"][i])
            self.assertEqual(dct_wave_function, wave_functions[i])
        self.assertEqual(len(dct["energy_levels"]), len(dct["wave_functions"]))
        dct_static_potential = spsc_data.Potential.from_dict(dct["static_potential"])
        self.assertEqual(dct_static_potential, static_potential)
        dct_density_potential = spsc_data.Potential.from_dict(dct["density_potential"])
        self.assertEqual(dct_density_potential, density_potential)
        dct_mass = spsc_data.MassValue.from_dict(dct["mass"])
        self.assertEqual(dct_mass, mass)

    def test_from_dict(self):
        energy_levels_number = random.randint(1, 5)
        rnd_energy_units = ElectronStateSimpleIOTestCase.random_units("EnergyValue")
        energy_levels = []
        for i in range(energy_levels_number):
            energy_levels.append(ElectronStateSimpleIOTestCase.random_value("EnergyValue", rnd_energy_units).to_dict())
        rnd_granularity = random.randint(10, 20)
        wave_functions = []
        for i in range(energy_levels_number):
            wave_functions.append(ElectronStateSimpleIOTestCase.random_value_array("WaveFunction", rnd_granularity).
                                  to_dict())
        static_potential = ElectronStateSimpleIOTestCase.random_value_array("Potential", rnd_granularity,
                                                                            rnd_energy_units).to_dict()
        density_potential = ElectronStateSimpleIOTestCase.random_value_array("Potential", rnd_granularity,
                                                                             rnd_energy_units).to_dict()
        mass = ElectronStateSimpleIOTestCase.random_value("MassValue").to_dict()
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

    # def test_import_export(self):
    #