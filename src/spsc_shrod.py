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

    @abstractmethod
    def _is_solution(self, solution):
        pass

    def solve(self, potential, mass, length):
        self.E_start.convert_to(self.dE.units)
        self.E_end.convert_to(self.dE.units)
        self.E_current = self.E_start
        self.iteration = self.iteration_factory.get_iteration(potential, mass, length)
        solutions = []
        while self.E_current.value < self.E_end.value:
            solution_candidate = self.iteration.solve(self.E_current)
            if self._is_solution(solution_candidate):
                solutions.append((self.E_current, solution_candidate[0]))
            self.solution_history.append(solution_candidate)

        return solutions


