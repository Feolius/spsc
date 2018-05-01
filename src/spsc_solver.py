import spsc_data
import spsc_shrod
from abc import ABCMeta, abstractproperty, abstractmethod
import numpy as np
import spsc_puass
import spsc_constants


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
        for i in xrange(len(solutions)):
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
        for i in xrange(len(solutions)):
            solution = solutions[i]
            self.state.electron_states[0].wave_functions.append(solution[1])
            self.state.electron_states[0].energy_levels.append(solution[0])


class LatticeSymmetrySolver(ASolver):
    def solve(self):
        # 30 nm well
        # E_start = spsc_data.EnergyValue(0.001, "eV")
        # E_end = spsc_data.EnergyValue(0.4, "eV")
        # 13 nm well
        E_start = spsc_data.EnergyValue(0.01, "eV")
        E_end = spsc_data.EnergyValue(0.12, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        static_potential = self.state.electron_states[0].static_potential
        meta_info = static_potential.meta_info
        mass = self.state.electron_states[0].mass
        if isinstance(mass, spsc_data.APhysValueArray):
            iteration_factory = spsc_shrod.SolutionIterationSymmetryLatticeDiffMassFactory()
        else:
            iteration_factory = spsc_shrod.SolutionIterationSymmetryLatticeFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategySymmetricWell(E_start, E_end, dE, 2, iteration_factory)
        density_potential = self.state.density_potential
        length = self.state.length
        electron_state = self.state.electron_states[0]
        for j in xrange(5):
            potential = static_potential + density_potential
            potential_offset = spsc_data.EnergyValue(potential[meta_info["well_start"]], potential.units)
            potential = potential - spsc_data.Potential(
                potential_offset.value * np.ones((len(potential),), "float64"), potential_offset.units)
            potential.meta_info = meta_info
            solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0))
            for i in xrange(len(solutions)):
                solution = solutions[i]
                if len(electron_state.wave_functions) > i:
                    electron_state.wave_functions[i] = solution[1]
                else:
                    electron_state.wave_functions.append(solution[1])
                energy_level = solution[0] + potential_offset
                if len(electron_state.energy_levels) > i:
                    electron_state.energy_levels[i] = energy_level
                else:
                    electron_state.energy_levels.append(energy_level)

            density_potential = self.__solve_puass()
        self.state.density_potential = density_potential

    def __get_E_Fermi(self):
        electron_state = self.state.electron_states[0]
        dos = self.__get_dos()
        # @TODO currently it's working for two and one subband systems only
        energy_levels_number = len(electron_state.energy_levels)
        if energy_levels_number > 2:
            raise StandardError("Cannot find E Fermi for more than 2 energy levels filled")
        for i in xrange(energy_levels_number):
            electron_state.energy_levels[i].convert_to(electron_state.energy_levels[i].units_default)

        electron_state.sum_density.convert_to(electron_state.sum_density.units_default)
        n = electron_state.sum_density.value
        if energy_levels_number == 1:
            E_Fermi_value = n / dos + electron_state.energy_levels[0].value
        else:
            E1 = electron_state.energy_levels[0].value
            E2 = electron_state.energy_levels[1].value
            E_Fermi_value = 0.5 * (n / dos + E1 + E2)
            if E_Fermi_value < E2:
                E_Fermi_value = n / dos + electron_state.energy_levels[0].value

        E_Fermi = spsc_data.EnergyValue(E_Fermi_value)
        # E_Fermi.convert_to("eV")
        return E_Fermi

    def __solve_puass(self):
        electron_state = self.state.electron_states[0]
        dos = self.__get_dos()
        # @TODO currently it's working for two and one subband systems only
        energy_levels_number = len(electron_state.energy_levels)
        if energy_levels_number > 2:
            raise StandardError("Cannot find E Fermi for more than 2 energy levels filled")
        E_Fermi = self.__get_E_Fermi()
        for i in xrange(energy_levels_number):
            electron_state.energy_levels[i].convert_to(electron_state.energy_levels[i].units_default)
        E_Fermi.convert_to(E_Fermi.units_default)
        for i in xrange(energy_levels_number):
            electron_state.energy_levels[i].convert_to(electron_state.energy_levels[i].units_default)
        h = 1.0 / (self.state.get_granularity() - 1)
        electron_state.sum_density.convert_to(electron_state.sum_density.units_default)
        if energy_levels_number == 1 or E_Fermi.value < electron_state.energy_levels[1].value:
            n = electron_state.sum_density.value
            electron_density = spsc_data.Density(
                n * h * (electron_state.wave_functions[0].value ** 2))
        else:
            n1 = dos * (E_Fermi.value - electron_state.energy_levels[0].value)
            n2 = dos * (E_Fermi.value - electron_state.energy_levels[1].value)
            density1 = spsc_data.Density(n1 * h * electron_state.wave_functions[0].value ** 2)
            density2 = spsc_data.Density(n2 * h * electron_state.wave_functions[1].value ** 2)
            electron_density = density1 + density2

        density = self.state.static_density - electron_density
        puass_solution_strategy = spsc_puass.GaussSolutionStrategy()
        return puass_solution_strategy.solve(density, 12, self.state.length)

    def __get_dos(self):
        self.state.electron_states[0].mass.convert_to(self.state.electron_states[0].mass.units_default)
        if isinstance(self.state.electron_states[0].mass, spsc_data.APhysValueArray):
            mass = spsc_data.MassValue(self.state.electron_states[0].mass[self.state.get_granularity() / 2])
        else:
            mass = self.state.electron_states[0].mass
        return mass.value / (np.pi * spsc_constants.h_plank ** 2)


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
        for j in xrange(10):
            potential = static_potential + density_potential
            potential_offset = spsc_data.EnergyValue(potential[meta_info["well_start"]], potential.units)
            potential = potential - spsc_data.Potential(
                potential_offset.value * np.ones((len(potential),), "float64"), potential_offset.units)
            potential.meta_info = meta_info
            solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0, 10.0 ** -25, 0))
            for i in xrange(len(solutions)):
                solution = solutions[i]
                if len(electron_state.wave_functions) > i:
                    electron_state.wave_functions[i] = solution[1]
                else:
                    electron_state.wave_functions.append(solution[1])
                energy_level = solution[0] + potential_offset
                if len(electron_state.energy_levels) > i:
                    electron_state.energy_levels[i] = energy_level
                else:
                    electron_state.energy_levels.append(energy_level)
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
        symmetry_solver = LatticeSymmetrySolver(self.state)
        symmetry_solver.solve()

        for j in xrange(5):
            self.__gamma_shrod()
            self.__x_shrod()
            self.state.density_potential = self.__solve_puass()

        # E_start = spsc_data.EnergyValue(0.005, "eV")
        # E_end = spsc_data.EnergyValue(0.3, "eV")
        # dE = spsc_data.EnergyValue(0.0001, "eV")
        #
        # meta_info = self.state.electron_states[1].static_potential.meta_info
        # self.state.electron_states[1].static_potential.convert_to('eV')
        # self.state.density_potential.convert_to('eV')
        # static_potential_arr = self.state.electron_states[1].static_potential.value[
        #                        meta_info['x_solution_start']:meta_info['x_solution_end']]
        # density_potential_arr = self.state.density_potential.value[
        #                         meta_info['x_solution_start']:meta_info['x_solution_end']]
        #
        # potential_arr = static_potential_arr + density_potential_arr
        # potential = spsc_data.Potential(potential_arr, "eV")
        # potential_offset = spsc_data.EnergyValue(np.amin(potential_arr), 'eV')
        # potential = potential - spsc_data.Potential(
        #     potential_offset.value * np.ones((len(potential),), "float64"), potential_offset.units)
        # potential.meta_info = meta_info
        #
        # mass = self.state.electron_states[1].mass
        # length = self.state.length
        # iteration_factory = spsc_shrod.SolutionIterationRungeKuttFactory()
        # solution_strategy = spsc_shrod.IterableSolutionStrategyNonSymmetricWell(E_start, E_end, dE, 1,
        #                                                                         iteration_factory)
        # solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0, 10.0 ** -25, 0))
        #
        # wave_function = solutions[0][1]
        # zeros = np.zeros((len(self.state.density_potential),))
        # wave_function_full_arr = np.concatenate((zeros[:meta_info['x_solution_start']], wave_function.value, zeros[meta_info['x_solution_end']:]))
        # wave_function = spsc_data.WaveFunction(wave_function_full_arr)
        # wave_function.mirror()
        # wave_function.normalize()
        # self.state.electron_states[1].wave_functions.append(wave_function)
        #
        # energy_level = solutions[0][0]
        # energy_level = energy_level + potential_offset
        # self.state.electron_states[1].energy_levels.append(energy_level)

    def __gamma_shrod(self):
        E_start = spsc_data.EnergyValue(0.001, "eV")
        E_end = spsc_data.EnergyValue(0.4, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")
        static_potential = self.state.electron_states[0].static_potential
        meta_info = static_potential.meta_info
        mass = self.state.electron_states[0].mass
        if isinstance(mass, spsc_data.APhysValueArray):
            iteration_factory = spsc_shrod.SolutionIterationSymmetryLatticeDiffMassFactory()
        else:
            iteration_factory = spsc_shrod.SolutionIterationSymmetryLatticeFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategySymmetricWell(E_start, E_end, dE, 2, iteration_factory)
        density_potential = self.state.density_potential
        length = self.state.length
        electron_state = self.state.electron_states[0]
        potential = static_potential + density_potential
        potential_offset = spsc_data.EnergyValue(potential[meta_info["well_start"]], potential.units)
        potential = potential - spsc_data.Potential(
            potential_offset.value * np.ones((len(potential),), "float64"), potential_offset.units)
        potential.meta_info = meta_info
        solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0))
        for i in xrange(len(solutions)):
            solution = solutions[i]
            if len(electron_state.wave_functions) > i:
                electron_state.wave_functions[i] = solution[1]
            else:
                electron_state.wave_functions.append(solution[1])
            energy_level = solution[0] + potential_offset
            if len(electron_state.energy_levels) > i:
                electron_state.energy_levels[i] = energy_level
            else:
                electron_state.energy_levels.append(energy_level)

    def __x_shrod(self):
        E_start = spsc_data.EnergyValue(0.005, "eV")
        E_end = spsc_data.EnergyValue(0.3, "eV")
        dE = spsc_data.EnergyValue(0.0001, "eV")

        meta_info = self.state.electron_states[1].static_potential.meta_info
        self.state.electron_states[1].static_potential.convert_to('eV')
        self.state.density_potential.convert_to('eV')
        static_potential_arr = self.state.electron_states[1].static_potential.value[
                               meta_info['x_solution_start']:meta_info['x_solution_end']]
        density_potential_arr = self.state.density_potential.value[
                                meta_info['x_solution_start']:meta_info['x_solution_end']]

        potential_arr = static_potential_arr + density_potential_arr
        potential = spsc_data.Potential(potential_arr, "eV")
        potential_offset = spsc_data.EnergyValue(np.amin(potential_arr), 'eV')
        potential = potential - spsc_data.Potential(
            potential_offset.value * np.ones((len(potential),), "float64"), potential_offset.units)
        potential.meta_info = meta_info

        mass = self.state.electron_states[1].mass
        length = self.state.length
        iteration_factory = spsc_shrod.SolutionIterationRungeKuttFactory()
        solution_strategy = spsc_shrod.IterableSolutionStrategyNonSymmetricWell(E_start, E_end, dE, 1,
                                                                                iteration_factory)
        solutions = solution_strategy.solve(potential, mass, length, (10.0 ** -20, 0, 10.0 ** -25, 0))

        wave_function = solutions[0][1]
        zeros = np.zeros((len(self.state.density_potential),))
        wave_function_full_arr = np.concatenate(
            (zeros[:meta_info['x_solution_start']], wave_function.value, zeros[meta_info['x_solution_end']:]))
        wave_function = spsc_data.WaveFunction(wave_function_full_arr)
        wave_function.mirror()
        wave_function.normalize()
        if len(self.state.electron_states[1].wave_functions) > 0:
            self.state.electron_states[1].wave_functions[0] = wave_function
        else:
            self.state.electron_states[1].wave_functions.append(wave_function)

        energy_level = solutions[0][0]
        energy_level = energy_level + potential_offset
        if len(self.state.electron_states[1].energy_levels) > 0:
            self.state.electron_states[1].energy_levels[0] = energy_level
        else:
            self.state.electron_states[1].energy_levels.append(energy_level)

    def __get_E_Fermi(self):
        electron_state_gamma = self.state.electron_states[0]
        electron_state_x = self.state.electron_states[1]
        dos_gamma = self.__get_dos_gamma()
        dos_x = self.__get_dos_x()
        # @TODO currently it's working for two and one subband systems only
        energy_levels_number = len(electron_state_gamma.energy_levels)
        if energy_levels_number > 2:
            raise StandardError("Cannot find E Fermi for more than 2 energy levels filled")
        for i in xrange(energy_levels_number):
            electron_state_gamma.energy_levels[i].convert_to(electron_state_gamma.energy_levels[i].units_default)
        electron_state_x.energy_levels[0].convert_to(electron_state_x.energy_levels[0].units_default)

        electron_state_gamma.sum_density.convert_to(electron_state_gamma.sum_density.units_default)
        n = electron_state_gamma.sum_density.value
        Ex = electron_state_x.energy_levels[0].value
        E1 = electron_state_gamma.energy_levels[0].value
        if energy_levels_number == 1:
            E_Fermi_value = (n + dos_gamma * E1 + dos_x * Ex) / (dos_gamma + dos_x)
            E_max = max(Ex, E1)
            if E_Fermi_value < E_max:
                if E_max == Ex:
                    dos = dos_gamma
                else:
                    dos = dos_x
                E_Fermi_value = n / dos + min(Ex, E1)
        else:
            E2 = electron_state_gamma.energy_levels[1].value
            E_Fermi_value = (n + dos_gamma * (E1 + E2) + dos_x * Ex) / (2 * dos_gamma + dos_x)
            E_max = max(Ex, E1, E2)
            if E_Fermi_value < E_max:
                if E_max == Ex:
                    E_Fermi_value = (n + dos_gamma * (E1 + E2)) / (2 * dos_gamma)
                    if E_Fermi_value < E2:
                        E_Fermi_value = n / dos_gamma + E1
                else:
                    E_Fermi_value = (n + dos_gamma * E1 + dos_x * Ex) / (dos_gamma + dos_x)
                    E_max = max(Ex, E1)
                    if E_Fermi_value < E_max:
                        if E_max == Ex:
                            dos = dos_gamma
                        else:
                            dos = dos_x
                        E_Fermi_value = n / dos + min(Ex, E1)
        E_Fermi = spsc_data.EnergyValue(E_Fermi_value)
        return E_Fermi

    def __solve_puass(self):
        electron_state_gamma = self.state.electron_states[0]
        electron_state_x = self.state.electron_states[1]
        dos_gamma = self.__get_dos_gamma()
        dos_x = self.__get_dos_x()
        # @TODO currently it's working for two and one subband systems only
        energy_levels_number = len(electron_state_gamma.energy_levels)
        if energy_levels_number > 2:
            raise StandardError("Cannot find E Fermi for more than 2 energy levels filled")
        E_Fermi = self.__get_E_Fermi()
        E_Fermi.convert_to(E_Fermi.units_default)
        for i in xrange(energy_levels_number):
            electron_state_gamma.energy_levels[i].convert_to(electron_state_gamma.energy_levels[i].units_default)
        electron_state_x.energy_levels[0].convert_to(electron_state_x.energy_levels[0].units_default)

        n1 = dos_gamma * (E_Fermi.value - electron_state_gamma.energy_levels[0].value)
        n2 = n_x = 0
        if energy_levels_number == 2 and E_Fermi.value > electron_state_gamma.energy_levels[1].value:
            n2 = dos_gamma * (E_Fermi.value - electron_state_gamma.energy_levels[1].value)
        if E_Fermi.value > electron_state_x.energy_levels[0].value:
            n_x = dos_x * (E_Fermi.value - electron_state_x.energy_levels[0].value)

        h = 1.0 / (self.state.get_granularity() - 1)
        electron_state_gamma.sum_density.convert_to(electron_state_gamma.sum_density.units_default)
        density1 = spsc_data.Density(n1 * h * electron_state_gamma.wave_functions[0].value ** 2)
        density2 = spsc_data.Density(n2 * h * electron_state_gamma.wave_functions[1].value ** 2)
        density_x = spsc_data.Density(n_x * h * electron_state_x.wave_functions[0].value ** 2)
        electron_density = density1 + density2 + density_x
        density = self.state.static_density - electron_density
        puass_solution_strategy = spsc_puass.GaussSolutionStrategy()
        return puass_solution_strategy.solve(density, 12, self.state.length)

    def __get_dos_gamma(self):
        self.state.electron_states[0].mass.convert_to(self.state.electron_states[0].mass.units_default)
        if isinstance(self.state.electron_states[0].mass, spsc_data.APhysValueArray):
            mass = spsc_data.MassValue(self.state.electron_states[0].mass[self.state.get_granularity() / 2])
        else:
            mass = self.state.electron_states[0].mass
        return mass.value / (np.pi * spsc_constants.h_plank ** 2)

    def __get_dos_x(self):
        self.state.electron_states[1].mass.convert_to(self.state.electron_states[0].mass.units_default)
        if isinstance(self.state.electron_states[1].mass, spsc_data.APhysValueArray):
            mass = spsc_data.MassValue(self.state.electron_states[1].mass[self.state.get_granularity() / 2])
        else:
            mass = self.state.electron_states[1].mass
        return mass.value / (np.pi * spsc_constants.h_plank ** 2)
