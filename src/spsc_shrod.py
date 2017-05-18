import spsc_data
from abc import ABCMeta, abstractproperty, abstractmethod


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
        self.potential = potential
        self.mass = mass
        self.length = length

    @abstractmethod
    def solve(self, E):
        pass


class ASolutionIterationFactory(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_iteration(self, potential, mass, length):
        pass


class SolutionIterationSingleWell(ASolutionIteration):

    def solve(self, E):
        pass


class SolutionIterationSingleWellFactory(ASolutionIterationFactory):

    def get_iteration(self, potential, mass, length):
        return SolutionIterationSingleWell(potential, mass, length)

