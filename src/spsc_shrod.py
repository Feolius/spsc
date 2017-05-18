import spsc_data
from abc import ABCMeta, abstractproperty, abstractmethod
import numpy as np
import spsc_constants as constants
import math


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
        E_current = self.E_start
        iteration = self.iteration_factory.get_iteration(potential, mass, length)
        solutions = []
        while E_current.value < self.E_end.value:
            solution_candidate = iteration.solve(E_current)
            if self._is_solution(solution_candidate):
                solutions.append((E_current, solution_candidate[0]))
            self.solution_history.append(solution_candidate)
            E_current += self.dE
            self._count += 1
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
        self.potential = potential.convert_to(potential.units_default)
        self.mass = mass.convert_to(potential.units_default)
        self.length = length.convert_to(potential.units_default)

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
        N = len(self.potential)
        solution = (spsc_data.WaveFunction(np.zeros((N,))), spsc_data.WaveFunction(np.zeros((N,))))
        gamma = 2 * self.mass * (self.length ** 2) / (constants.h_plank ** 2)
        (A, B) = self._get_initial_AB(E, solution_start)
        potential_step = self.potential[0]
        for i in range(N):
            x = i / (N - 1)
            if self.potential[i] == potential_step:
              if E < potential_step:
                  k = math.sqrt((potential_step - E) * gamma)



        return solution

    def _get_initial_AB(self, E, solution_start):
        # TODO need to implement this method properly
        return solution_start[0], 0





class SolutionIterationFlatPotentialFactory(ASolutionIterationFactory):

    def get_iteration(self, potential, mass, length):
        return SolutionIterationFlatPotential(potential, mass, length)

