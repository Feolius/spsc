import spsc_data
import spsc_io
import spsc_shrod
from abc import ABCMeta, abstractproperty, abstractmethod


class ASolver(object):
    __metaclass__ = ABCMeta

    def __init__(self, state):
        self.state = state

    @abstractmethod
    def solve(self):
        pass


class ShrodSolverSimple(ASolver):

    def solve(self):
        E_start = spsc_data.EnergyValue(0.001, "eV")
        E_end = spsc_data.EnergyValue(1, "eV")
        dE = spsc_data.EnergyValue(0.001, "eV")
        iteration_factory = spsc_shrod.SolutionIterationFlatPotentialFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategySymmetricWell(E_start, E_end, dE, iteration_factory)
        potential = self.state.electron_states[0].static_potential
        mass = self.state.electron_states[0].mass
        length = self.state.length
        solution = solution_strategy.solve(potential, mass, length)



class AState(spsc_io.Default):
    __metaclass__ = ABCMeta

    def electron_states_getter(self):
        pass

    def electron_states_setter(self, electron_states):
        pass

    electron_states = abstractproperty(electron_states_getter, electron_states_setter)

    def static_density_getter(self):
        pass

    def static_density_setter(self, static_density):
        pass

    static_density = abstractproperty(static_density_getter, static_density_setter)

    def density_potential_getter(self):
        pass

    def density_potential_setter(self, density_potential):
        pass

    density_potential = abstractproperty(density_potential_getter, density_potential_setter)

    def length_getter(self):
        pass

    def length_setter(self, length):
        pass

    length = abstractproperty(length_getter, length_setter)


class StateSimple(AState):

    _electron_states = []

    def electron_states_getter(self):
        return self._electron_states

    def electron_states_setter(self, electron_states):
        self._electron_states = electron_states

    electron_states = property(electron_states_getter, electron_states_setter)

    _static_density = spsc_data.Density([])

    def static_density_getter(self):
        return self._static_density

    def static_density_setter(self, static_density):
        self._static_density = static_density

    static_density = property(static_density_getter, static_density_setter)

    _density_potential = spsc_data.Potential([])

    def density_potential_getter(self):
        return self._density_potential

    def density_potential_setter(self, density_potential):
        self._density_potential = density_potential

    density_potential = property(density_potential_getter, density_potential_setter)

    _length = spsc_data.LengthValue(0)

    def length_getter(self):
        return self._length

    def length_setter(self, length):
        self._length = length

    length = property(length_getter, length_setter)

    @classmethod
    def from_dict(cls, dct):
        state = cls()
        if "electron_states" in dct:
            electron_states = []
            for electron_state_dct in dct["electron_states"]:
                electron_states.append(ElectronStateSimple.from_dict(electron_state_dct))
            state.electron_states = electron_states
        if "static_density" in dct:
            state.static_density = spsc_data.Density.from_dict(dct["static_density"])
        if "density_potential" in dct:
            state.density_potential = spsc_data.Potential.from_dict(dct["density_potential"])
        if "length" in dct:
            state.length = spsc_data.LengthValue.from_dict(dct["length"])
        if not state._is_granularity_consistent():
            raise ImportError("Electron State data have different granularity in State dump")
        return state

    def to_dict(self):
        dct = {
            "electron_states": [],
            "static_density": self.static_density.to_dict(),
            "density_potential": self.density_potential.to_dict(),
            "length": self.length.to_dict()
        }
        for electron_state in self.electron_states:
            dct["electron_states"].append(electron_state.to_dict())
        return dct

    def _is_granularity_consistent(self):
        consistent = True
        if len(self.electron_states) > 0:
            granularity = self.electron_states[0].get_granularity()
            for i in range(len(self.electron_states) - 1):
                consistent = consistent and granularity == self.electron_states[i + 1].get_granularity()
        return consistent

    def get_granularity(self):
        granularity = 0
        if self._is_granularity_consistent():
            if self.electron_states:
                granularity = self.electron_states[0].get_granularity()
        else:
            raise StandardError("State granularity is inconsistent")
        return granularity

    def __eq__(self, other):
        equal = False
        if type(other) == self.__class__:
            electron_states_equal = False
            if len(self.electron_states) == len(other.electron_states):
                electron_states_equal = True
                for i in range(len(self.electron_states)):
                    electron_states_equal = electron_states_equal and self.electron_states[i] == other.electron_states[i]
            static_density_equal = self.static_density == other.static_density
            density_potential_equal = self.density_potential == other.density_potential
            length_equal = self.length == other.length
            equal = electron_states_equal and static_density_equal and density_potential_equal and length_equal
        return equal

    def __ne__(self, other):
        return not self == other


class AElectronState(spsc_io.Default):
    __metaclass__ = ABCMeta

    def wave_functions_getter(self):
        pass

    def wave_functions_setter(self, wave_functions):
        pass

    wave_functions = abstractproperty(wave_functions_getter, wave_functions_setter)

    def energy_levels_getter(self):
        pass

    def energy_levels_setter(self, energy_levels):
        pass

    energy_levels = abstractproperty(energy_levels_getter, energy_levels_setter)

    def static_potential_getter(self):
        pass

    def static_potential_setter(self, static_potential):
        pass

    static_potential = abstractproperty(static_potential_getter, static_potential_setter)

    def mass_getter(self):
        pass

    def mass_setter(self, mass):
        pass

    mass = abstractproperty(mass_getter, mass_setter)


class ElectronStateSimple(AElectronState):
    _wave_functions = []

    def wave_functions_getter(self):
        return self._wave_functions

    def wave_functions_setter(self, wave_functions):
        self._wave_functions = wave_functions

    wave_functions = property(wave_functions_getter, wave_functions_setter)

    _energy_levels = []

    def energy_levels_getter(self):
        return self._energy_levels

    def energy_levels_setter(self, energy_levels):
        self._energy_levels = energy_levels

    energy_levels = property(energy_levels_getter, energy_levels_setter)

    _static_potential = spsc_data.Potential([])

    def static_potential_getter(self):
        return self._static_potential

    def static_potential_setter(self, static_potential):
        self._static_potential = static_potential

    static_potential = property(static_potential_getter, static_potential_setter)

    _density_potential = spsc_data.Potential([])

    _mass = spsc_data.MassValue(0)

    def mass_getter(self):
        return self._mass

    def mass_setter(self, mass):
        self._mass = mass

    mass = property(mass_getter, mass_setter)

    @classmethod
    def from_dict(cls, dct):
        electron_sate = cls()
        if "static_potential" in dct:
            electron_sate.static_potential = spsc_data.Potential.from_dict(dct["static_potential"])
        if "wave_functions" in dct and type(dct["wave_functions"]) is list:
            wave_functions = []
            for wave_function_dct in dct["wave_functions"]:
                wave_function = spsc_data.WaveFunction.from_dict(wave_function_dct)
                wave_functions.append(wave_function)
            electron_sate.wave_functions = wave_functions
        if "energy_levels" in dct and type(dct["energy_levels"]) is list:
            energy_levels = []
            for energy_level in dct["energy_levels"]:
                energy_levels.append(spsc_data.EnergyValue.from_dict(energy_level))
            electron_sate.energy_levels = energy_levels
        if "mass" in dct:
            electron_sate.mass = spsc_data.MassValue.from_dict(dct["mass"])
        if not electron_sate._is_granularity_consistent():
            raise ImportError("Incoming data arrays have different lengths in Electron State dump.")
        if not electron_sate._is_wave_functions_consistent():
            raise ImportError("Number of energy levels and wave functions doesn't match in Electron State dump.")
        return electron_sate

    def to_dict(self):
        dct = {
            "static_potential": self.static_potential.to_dict(),
            "wave_functions": [],
            "energy_levels": [],
            "mass": self.mass.to_dict()
        }
        for wave_function in self.wave_functions:
            dct["wave_functions"].append(wave_function.to_dict())
        for energy_level in self.energy_levels:
            dct["energy_levels"].append(energy_level.to_dict())
        return dct

    def __eq__(self, other):
        equal = False
        if type(other) == self.__class__:
            wave_functions_equal = False
            if len(self.wave_functions) == len(other.wave_functions):
                wave_functions_equal = True
                for i in range(len(self.wave_functions)):
                    wave_functions_equal = wave_functions_equal and self.wave_functions[i] == other.wave_functions[i]
            energy_levels_equal = False
            if len(self.energy_levels) == len(other.energy_levels):
                energy_levels_equal = True
                for i in range(len(self.energy_levels)):
                    energy_levels_equal = energy_levels_equal and self.energy_levels[i] == other.energy_levels[i]
            static_potential_equal = self.static_potential == other.static_potential
            mass_equal = self.mass == other.mass
            equal = wave_functions_equal and energy_levels_equal and static_potential_equal and mass_equal
        return equal

    def __ne__(self, other):
        return not self == other

    def _is_granularity_consistent(self):
        granularity = len(self.static_potential)
        consistent = True
        for wave_function in self.wave_functions:
            consistent = consistent and granularity == len(wave_function)
        return consistent

    def _is_wave_functions_consistent(self):
        return len(self.wave_functions) == len(self.energy_levels)

    def get_granularity(self):
        if self._is_granularity_consistent():
            granularity = len(self.static_potential)
        else:
            raise StandardError("Electron state granularity is inconsistent")
        return granularity


