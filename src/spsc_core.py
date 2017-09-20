import spsc_data
import spsc_io
import spsc_shrod
from abc import ABCMeta, abstractproperty, abstractmethod
import matplotlib.pyplot as plt
import numpy as np
import spsc_puass


class ASolver(object):
    __metaclass__ = ABCMeta

    def __init__(self, state):
        self.state = state

    @abstractmethod
    def solve(self):
        pass


class SymmetricWellSolver(ASolver):
    def solve(self):
        E_start = spsc_data.EnergyValue(0.001, "eV")
        E_end = spsc_data.EnergyValue(0.4, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        iteration_factory = spsc_shrod.SolutionIterationFlatPotentialFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategySymmetricWell(E_start, E_end, dE, 6, iteration_factory)
        potential = self.state.electron_states[0].static_potential
        mass = self.state.electron_states[0].mass
        length = self.state.length
        solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0))
        for i in range(len(solutions)):
            solution = solutions[i]
            self.state.electron_states[0].wave_functions.append(solution[1])
            self.state.electron_states[0].energy_levels.append(solution[0])


class SlopedWellSolver(ASolver):
    def solve(self):
        E_start = spsc_data.EnergyValue(0.03, "eV")
        E_end = spsc_data.EnergyValue(0.6, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        iteration_factory = spsc_shrod.SolutionIterationSlopePotentialFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategyNonSymmetricWell(E_start, E_end, dE, 8,
                                                                                iteration_factory)
        potential = self.state.electron_states[0].static_potential
        mass = self.state.electron_states[0].mass
        length = self.state.length
        solutions = solution_strategy.solve(potential, mass, length,
                                            (10.0 ** -20, 10.0 ** -15, 10.0 ** -25, -10.0 ** -15))
        for i in range(len(solutions)):
            solution = solutions[i]
            self.state.electron_states[0].wave_functions.append(solution[1])
            self.state.electron_states[0].energy_levels.append(solution[0])


class LatticeSymmetrySolver(ASolver):
    def solve(self):
        E_start = spsc_data.EnergyValue(0.001, "eV")
        E_end = spsc_data.EnergyValue(0.4, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        static_potential = self.state.electron_states[0].static_potential
        meta_info = static_potential.meta_info
        iteration_factory = spsc_shrod.SolutionIterationSymmetryLatticeFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategySymmetricWell(E_start, E_end, dE, 1, iteration_factory)
        density_potential = self.state.density_potential
        mass = self.state.electron_states[0].mass
        length = self.state.length
        electron_state = self.state.electron_states[0]
        for j in range(10):
            potential = static_potential + density_potential
            potential = potential - spsc_data.Potential(
                potential[meta_info["well_start"]] * np.ones((len(potential),), "float64"), potential.units)
            potential.meta_info = meta_info
            solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0))
            for i in range(len(solutions)):
                solution = solutions[i]
                if len(electron_state.wave_functions) > i:
                    electron_state.wave_functions[i] = solution[1]
                else:
                    electron_state.wave_functions.append(solution[1])
                if len(electron_state.energy_levels) > i:
                    electron_state.energy_levels[i] = solution[0]
                else:
                    electron_state.energy_levels.append(solution[0])
            h = 1.0 / (len(self.state.electron_states[0].wave_functions[0]) - 1)
            electron_density = spsc_data.Density(
                electron_state.sum_density.value * h * (electron_state.wave_functions[0].value ** 2),
                electron_state.sum_density.units)
            self.state.static_density.convert_to("m^-2")
            density = self.state.static_density - electron_density
            puass_solution_strategy = spsc_puass.GaussSolutionStrategy()
            density_potential = puass_solution_strategy.solve(density, 12, self.state.length)
        self.state.density_potential = density_potential


class LatticeSlopedSolver(ASolver):
    def solve(self):
        E_start = spsc_data.EnergyValue(0.0001, "eV")
        E_end = spsc_data.EnergyValue(0.2, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        static_potential = self.state.electron_states[0].static_potential
        meta_info = static_potential.meta_info
        iteration_factory = spsc_shrod.SolutionIterationSlopedLatticeFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategyNonSymmetricWell(E_start, E_end, dE, 1,
                                                                                iteration_factory)
        density_potential = self.state.density_potential
        mass = self.state.electron_states[0].mass
        length = self.state.length
        electron_state = self.state.electron_states[0]
        for j in range(10):
            potential = static_potential + density_potential
            potential = potential - spsc_data.Potential(
                potential[meta_info["well_start"]] * np.ones((len(potential),), "float64"), potential.units)
            potential.meta_info = meta_info
            solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0, 10.0 ** -25, 0))
            for i in range(len(solutions)):
                solution = solutions[i]
                if len(electron_state.wave_functions) > i:
                    electron_state.wave_functions[i] = solution[1]
                else:
                    electron_state.wave_functions.append(solution[1])
                if len(electron_state.energy_levels) > i:
                    electron_state.energy_levels[i] = solution[0]
                else:
                    electron_state.energy_levels.append(solution[0])
            h = 1.0 / (len(self.state.electron_states[0].wave_functions[0]) - 1)
            electron_density = spsc_data.Density(
                electron_state.sum_density.value * h * (electron_state.wave_functions[0].value ** 2),
                electron_state.sum_density.units)
            self.state.static_density.convert_to("m^-2")
            density = self.state.static_density - electron_density
            puass_solution_strategy = spsc_puass.GaussSolutionStrategy()
            prev_density_potential = density_potential
            density_potential = puass_solution_strategy.solve(density, 12, self.state.length)
            density_potential.convert_to(prev_density_potential.units)
            density_potential.value = (density_potential.value + prev_density_potential.value) / 2
        self.state.density_potential = density_potential


class LatticeXElectronsSymmetrySolver(ASolver):
    def solve(self):
        # symmetry_solver = LatticeSymmetrySolver(self.state)
        # symmetry_solver.solve()

        E_start = spsc_data.EnergyValue(0.01, "eV")
        E_end = spsc_data.EnergyValue(0.2, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")

        meta_info = self.state.electron_states[1].static_potential.meta_info
        self.state.electron_states[1].static_potential.convert_to('eV')
        self.state.density_potential.convert_to('eV')
        static_potential_arr = self.state.electron_states[1].static_potential.value[
                               meta_info['x_solution_start']:meta_info['x_solution_end']]
        density_potential_arr = self.state.density_potential.value[
                                meta_info['x_solution_start']:meta_info['x_solution_end']]

        potential_arr = static_potential_arr + density_potential_arr
        potential_arr = potential_arr - np.amin(potential_arr)
        potential = spsc_data.Potential(potential_arr, "eV")

        mass = self.state.electron_states[1].mass
        length = self.state.length
        iteration_factory = spsc_shrod.SolutionIterationRungeKuttFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategyNonSymmetricWell(E_start, E_end, dE, 1,
                                                                                iteration_factory)
        solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0, 10.0 ** -25, 0))
        wave_function = solutions[0][1]
        wave_function.instant_plot()


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
    _electron_states = None

    def electron_states_getter(self):
        return self._electron_states

    def electron_states_setter(self, electron_states):
        self._electron_states = electron_states

    electron_states = property(electron_states_getter, electron_states_setter)

    _static_density = None

    def static_density_getter(self):
        return self._static_density

    def static_density_setter(self, static_density):
        self._static_density = static_density

    static_density = property(static_density_getter, static_density_setter)

    _density_potential = None

    def density_potential_getter(self):
        return self._density_potential

    def density_potential_setter(self, density_potential):
        self._density_potential = density_potential

    density_potential = property(density_potential_getter, density_potential_setter)

    _length = None

    def length_getter(self):
        return self._length

    def length_setter(self, length):
        self._length = length

    length = property(length_getter, length_setter)

    def __init__(self):
        self.electron_states = []
        self.static_density = spsc_data.Density([])
        self.density_potential = spsc_data.Potential([])
        self.length = spsc_data.LengthValue(0)

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
                    electron_states_equal = electron_states_equal and self.electron_states[i] == other.electron_states[
                        i]
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

    def static_sum_density_getter(self):
        pass

    def static_sum_density_setter(self, sum_density):
        pass

    sum_density = abstractproperty(static_sum_density_getter, static_sum_density_setter)


class ElectronStateSimple(AElectronState):
    _wave_functions = None

    def wave_functions_getter(self):
        return self._wave_functions

    def wave_functions_setter(self, wave_functions):
        self._wave_functions = wave_functions

    wave_functions = property(wave_functions_getter, wave_functions_setter)

    _energy_levels = None

    def energy_levels_getter(self):
        return self._energy_levels

    def energy_levels_setter(self, energy_levels):
        self._energy_levels = energy_levels

    energy_levels = property(energy_levels_getter, energy_levels_setter)

    _static_potential = None

    def static_potential_getter(self):
        return self._static_potential

    def static_potential_setter(self, static_potential):
        self._static_potential = static_potential

    static_potential = property(static_potential_getter, static_potential_setter)

    _mass = None

    def mass_getter(self):
        return self._mass

    def mass_setter(self, mass):
        self._mass = mass

    mass = property(mass_getter, mass_setter)

    _sum_density = None

    def sum_density_getter(self):
        return self._sum_density

    def sum_density_setter(self, sum_density):
        self._sum_density = sum_density

    sum_density = property(sum_density_getter, sum_density_setter)

    def __init__(self):
        self.wave_functions = []
        self.energy_levels = []
        self.static_potential = spsc_data.Potential([])
        self.mass = spsc_data.MassValue(0)
        self.sum_density = spsc_data.DensityValue(0)

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
        if "sum_density" in dct:
            electron_sate.sum_density = spsc_data.DensityValue.from_dict(dct["sum_density"])
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
            "mass": self.mass.to_dict(),
            "sum_density": self.sum_density.to_dict()
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
