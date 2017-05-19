import spsc_data
from abc import ABCMeta, abstractmethod
import numpy as np
import spsc_constants as constants
import copy
import matplotlib.pyplot as plt


class ASolutionStrategy(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def solve(self, potential, mass, length):
        pass


class AIterableSolutionStrategy(ASolutionStrategy):
    __metaclass__ = ABCMeta

    def __init__(self, E_start, E_end, dE, iteration_factory):
        self.E_start = E_start
        self.E_end = E_end
        self.dE = dE
        self.iteration_factory = iteration_factory
        self.solution_history = []
        self._count = 0

    @abstractmethod
    def _is_solution(self, solution):
        pass

    def solve(self, potential, mass, length):
        self.E_start.convert_to(self.dE.units)
        self.E_end.convert_to(self.dE.units)
        E_current = copy.deepcopy(self.E_start)
        iteration = self.iteration_factory.get_iteration(potential, mass, length)
        solutions = []
        while E_current.value < self.E_end.value:
            solution_candidate = iteration.solve(E_current, (10.0 ** -20, 0))
            if self._is_solution(solution_candidate):
                solutions.append((E_current, solution_candidate[0]))
            self.solution_history.append(solution_candidate)
            plt.gcf().clear()
            plt.plot(solution_candidate[0].value)
            plt.ion()
            plt.pause(1)
            plt.show()
            solution_candidate[0].instant_plot()
            E_current += self.dE
            self._count += 1
            print self._count
        return solutions


class IterableSolutionStrategySymmetricWell(AIterableSolutionStrategy):
    def _is_solution(self, solution):
        is_solution = False
        if self.solution_history:
            prev_solution = self.solution_history[-1]
            middle_index = len(solution[0]) / 2
            middle_der = solution[1][middle_index]
            prev_middle_der = prev_solution[1][middle_index]
            is_der_changed = False
            if middle_der * prev_middle_der < 0:
                is_der_changed = True
            middle_func = solution[0][middle_index]
            prev_middle_func = prev_solution[0][middle_index]
            is_func_changed = False
            if middle_func * prev_middle_func < 0:
                is_func_changed = True
            if is_func_changed and is_der_changed:
                raise StandardError("Both function and derivation is changed in current iteration.")
            if is_func_changed or is_der_changed:
                is_solution = True
        return is_solution


class ASolutionIteration(object):
    __metaclass__ = ABCMeta

    def __init__(self, potential, mass, length):
        potential.convert_to(potential.units_default)
        mass.convert_to(mass.units_default)
        length.convert_to(length.units_default)
        self.potential = potential.value
        self.mass = mass.value
        self.length = length.value

    @abstractmethod
    def solve(self, E, solution_start):
        pass


class ASolutionIterationFactory(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_iteration(self, potential, mass, length):
        pass


class SolutionIterationFlatPotential(ASolutionIteration):
    def solve(self, E, solution_start):
        E.convert_to(E.units_default)
        E = E.value
        N = len(self.potential)
        solution = (spsc_data.WaveFunction(np.zeros((N,))), spsc_data.WaveFunction(np.zeros((N,))))
        gamma = 2.0 * self.mass * (self.length ** 2) / (constants.h_plank ** 2)
        (A, B) = self._get_initial_AB(E, solution_start)
        potential_level = self.potential[0]
        for i in range(N * 2 / 3):
            x = float(i) / (N - 1)
            if self.potential[i] != potential_level:
                AB_transition_matrix = self._get_AB_transition_matrix(potential_level, self.potential[i], E, x)
                (A, B) = np.dot(AB_transition_matrix, np.array([A, B]))
                potential_level = self.potential[i]
            if E < potential_level:
                k = np.sqrt((potential_level - E) * gamma)
                solution[0][i] = A * np.exp(k * x) + B * np.exp(-k * x)
                solution[1][i] = A * k * np.exp(k * x) - B * k * np.exp(-k * x)
            else:
                k = np.sqrt((E - potential_level) * gamma)
                solution[0][i] = A * np.sin(k * x) + B * np.cos(k * x)
                solution[1][i] = A * k * np.cos(k * x) - B * k * np.sin(k * x)

        return solution

    def _get_initial_AB(self, E, solution_start):
        # TODO need to implement this method properly
        return solution_start[0], 0

    def _get_AB_transition_matrix(self, prev_potential_level, new_potential_level, E, x):
        transition_matrix = np.zeros((2, 2))
        gamma = 2 * self.mass * (self.length ** 2) / (constants.h_plank ** 2)
        if E < prev_potential_level and E < new_potential_level:
            k1 = np.sqrt((prev_potential_level - E) * gamma)
            k2 = np.sqrt((new_potential_level - E) * gamma)
            transition_matrix[0][0] = np.exp(k1 * x) * np.exp(-k2 * x) * ((k1 / k2) + 1) / 2
            transition_matrix[0][1] = np.exp(-k1 * x) * np.exp(-k2 * x) * (1 - (k1 / k2)) / 2
            transition_matrix[1][0] = np.exp(k1 * x) * np.exp(k2 * x) * (1 - (k1 / k2)) / 2
            transition_matrix[1][1] = np.exp(-k1 * x) * np.exp(k2 * x) * ((k1 / k2) + 1) / 2
        elif E > prev_potential_level:
            k1 = np.sqrt((E - prev_potential_level) * gamma)
            k2 = np.sqrt((new_potential_level - E) * gamma)
            transition_matrix[0][0] = (k1 / (2 * k2) * np.cos(k1 * x) + np.sin(k1 * x) / 2) * np.exp(-k2 * x)
            transition_matrix[0][1] = (-k1 / (2 * k2) * np.sin(k1 * x) + np.cos(k1 * x) / 2) * np.exp(-k2 * x)
            transition_matrix[1][0] = (-k1 / (2 * k2) * np.cos(k1 * x) + np.sin(k1 * x) / 2) * np.exp(k2 * x)
            transition_matrix[1][1] = (k1 / (2 * k2) * np.sin(k1 * x) + np.cos(k1 * x) / 2) * np.exp(k2 * x)
        else:
            k1 = np.sqrt((prev_potential_level - E) * gamma)
            k2 = np.sqrt((E - new_potential_level) * gamma)
            transition_matrix[0][0] = np.exp(k1 * x) * (k1 / k2 * np.cos(k2 * x) + np.sin(k2 * x))
            transition_matrix[0][1] = np.exp(-k1 * x) * (np.sin(k2 * x) - k1 / k2 * np.cos(k2 * x))
            transition_matrix[1][0] = np.exp(k1 * x) * (np.cos(k2 * x) - k1 / (k2) * np.sin(k2 * x))
            transition_matrix[1][1] = np.exp(-k1 * x) * (k1 / k2 * np.sin(k2 * x) + np.cos(k2 * x))
        return transition_matrix


class SolutionIterationFlatPotentialFactory(ASolutionIterationFactory):
    def get_iteration(self, potential, mass, length):
        return SolutionIterationFlatPotential(potential, mass, length)
