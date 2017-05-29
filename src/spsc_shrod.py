import spsc_data
from abc import ABCMeta, abstractmethod
import numpy as np
import spsc_constants as constants
import copy
import scipy.special as special
import matplotlib.pyplot as plt


class ASolutionStrategy(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def solve(self, potential, mass, length, solution_start):
        pass


class AIterableSolutionStrategy(ASolutionStrategy):
    __metaclass__ = ABCMeta

    def __init__(self, E_start, E_end, dE, solutions_limit, iteration_factory):
        self.E_start = E_start
        self.E_end = E_end
        self.dE = dE
        self.solutions_limit = solutions_limit
        self.iteration_factory = iteration_factory
        self.solution_history = []
        self._count = 0
        self._w = []

    @abstractmethod
    def _is_solution(self, solution_candidate):
        pass

    @abstractmethod
    def _prepare_wave_function(self, solution_candidate):
        pass

    def solve(self, potential, mass, length, solution_start):
        self.E_start.convert_to(self.dE.units)
        self.E_end.convert_to(self.dE.units)
        E_current = copy.deepcopy(self.E_start)
        iteration = self.iteration_factory.get_iteration(potential, mass, length)
        solutions = []
        while E_current.value < self.E_end.value:
            solution_candidate = iteration.solve(E_current, solution_start)
            N = len(solution_candidate[0])
            # zeros = np.zeros((N,))
            # plt.gcf().clear()
            # a = np.concatenate((solution_candidate[0][0:(N - 1) * 2 / 3 + 1], zeros[(N - 1) * 2 / 3 + 1:]))
            # plt.ion()
            # plt.plot(a)
            # plt.show()
            # plt.pause(0.005)

            # plt.gcf().clear()
            # a = np.concatenate((zeros[:(N - 1) / 3 + 1], solution_candidate[2][(N - 1) / 3 + 1:]))
            # plt.ion()
            # plt.plot(a)
            # plt.show()
            # plt.pause(0.01)
            if self._is_solution(solution_candidate):
                wave_function = self._prepare_wave_function(solution_candidate)
                plt.gcf().clear()
                plt.plot(wave_function)
                potential.convert_to("eV")
                plt.plot(potential.value)
                E_current.convert_to("eV")
                print "Energy:", E_current.value
                plt.show()
                solutions.append((E_current, wave_function))

                if len(solutions) == self.solutions_limit:
                    break
            if self.solution_history and solution_candidate[0][N * 2 / 3] != self.solution_history[-1][0][N * 2 / 3]:
                solution_start = (solution_start[0], solution_start[1], -solution_start[2], solution_start[3])


            self.solution_history.append(solution_candidate)
            E_current.convert_to("eV")
            print E_current
            E_current.convert_to(self.dE.units)
            E_current += self.dE
            self._count += 1
        return solutions


class IterableSolutionStrategySymmetricWell(AIterableSolutionStrategy):
    def _is_solution(self, solution_candidate):
        is_solution = False
        if self.solution_history:
            prev_solution = self.solution_history[-1]
            middle_index = len(solution_candidate[0]) / 2
            middle_der = solution_candidate[1][middle_index]
            prev_middle_der = prev_solution[1][middle_index]
            is_der_changed = False
            if middle_der * prev_middle_der < 0:
                is_der_changed = True
            middle_func = solution_candidate[0][middle_index]
            prev_middle_func = prev_solution[0][middle_index]
            is_func_changed = False
            if middle_func * prev_middle_func < 0:
                is_func_changed = True
            if is_func_changed and is_der_changed:
                raise StandardError("Both function and derivation is changed in current iteration.")
            if is_func_changed or is_der_changed:
                is_solution = True
        return is_solution

    def _prepare_wave_function(self, solution_candidate):
        wave_function = solution_candidate[0].mirror()
        wave_function.normalize()
        return wave_function


class IterableSolutionStrategyNonSymmetricWell(AIterableSolutionStrategy):
    def _is_solution(self, solution_candidate):
        is_solution = False
        if self.solution_history:
            prev_solution = self.solution_history[-1]
            middle_index = len(solution_candidate[0]) / 2
            prev_w = prev_solution[1][middle_index] * prev_solution[2][middle_index] - \
                     prev_solution[3][middle_index] * prev_solution[0][middle_index]
            w = solution_candidate[1][middle_index] * solution_candidate[2][middle_index] - \
                solution_candidate[3][middle_index] * solution_candidate[0][middle_index]
            print "w:", w
            self._w.append(w)
            if w * prev_w < 0:
                is_solution = True
        return is_solution

    def _prepare_wave_function(self, solution_candidate):
        middle_index = len(solution_candidate[0]) / 2
        k = solution_candidate[0][middle_index] / solution_candidate[2][middle_index]
        wave_function_arr = np.concatenate(
            (solution_candidate[0].value[:middle_index + 1], k * solution_candidate[2].value[middle_index + 1:]))
        wave_function = spsc_data.WaveFunction(wave_function_arr)
        wave_function.normalize()
        return wave_function


class ASolutionIteration(object):
    __metaclass__ = ABCMeta

    def __init__(self, potential, mass, length):
        self.potential = potential
        self.mass = mass
        self.length = length

    @abstractmethod
    def solve(self, E, solution_start):
        pass

    def _reset_to_default_units(self):
        self.potential.convert_to(self.potential.units_default)
        self.mass.convert_to(self.mass.units_default)
        self.length.convert_to(self.length.units_default)


class ASolutionIterationFactory(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_iteration(self, potential, mass, length):
        pass


class SolutionIterationFlatPotential(ASolutionIteration):
    def solve(self, E, solution_start):
        self._reset_to_default_units()
        E.convert_to(E.units_default)
        E = E.value
        N = len(self.potential)
        solution = (spsc_data.WaveFunction(np.zeros((N,))), spsc_data.WaveFunction(np.zeros((N,))))
        gamma = 2.0 * self.mass.value * (self.length.value ** 2) / (constants.h_plank ** 2)
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
        gamma = 2 * self.mass.value * (self.length.value ** 2) / (constants.h_plank ** 2)
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


class SolutionIterationSlopePotential(ASolutionIteration):
    SLOPE_RESOLUTION_DOTS_NUM = 100

    def __init__(self, potential, mass, length):
        super(SolutionIterationSlopePotential, self).__init__(potential, mass, length)
        self.slope = spsc_data.ElectricFieldValue(0)
        self.x0 = spsc_data.LengthValue(0)
        self.eps0 = spsc_data.EnergyValue(0)
        self._reset_to_default_units()
        # Calculate slope
        energy_difference = self.potential[SolutionIterationSlopePotential.SLOPE_RESOLUTION_DOTS_NUM - 1] - \
                            self.potential[0]
        energy_difference = spsc_data.EnergyValue(energy_difference)
        voltage_difference = energy_difference.value / constants.e
        length_piece = self.length.value * (SolutionIterationSlopePotential.SLOPE_RESOLUTION_DOTS_NUM - 1) / (len(self.potential) - 1)
        self.slope.value = voltage_difference / length_piece

        self.x0.value = np.power(constants.h_plank ** 2 / (2 * self.mass.value * constants.e * self.slope.value), float(1) / 3)
        self.eps0.value = constants.e * self.slope.value * self.x0.value

    def _reset_to_default_units(self):
        super(SolutionIterationSlopePotential, self)._reset_to_default_units()
        self.slope.convert_to(self.slope.units_default)
        self.x0.convert_to(self.x0.units_default)
        self.eps0.convert_to(self.eps0.units_default)

    def solve(self, E, solution_start):
        self._reset_to_default_units()
        E.convert_to(E.units_default)
        N = len(self.potential)
        solution = (spsc_data.WaveFunction(np.zeros((N,))), spsc_data.WaveFunction(np.zeros((N,))),
                    spsc_data.WaveFunction(np.zeros((N,))), spsc_data.WaveFunction(np.zeros((N,))))
        potential_threshold = abs(2.0 * constants.e * self.slope.value * self.length.value / (len(self.potential) - 1))

        (A, B) = self._get_initial_left_AB(E, solution_start[0], solution_start[1])
        U = self._get_U(0)
        U.convert_to("eV")
        for i in range(N * 2 / 3):
            if i > 0 and abs(self.potential[i] - self.potential[i-1]) > potential_threshold:
                (A, B) = self._get_AB(E, i - 1, i, solution[0][i - 1], solution[1][i - 1])
                U = self._get_U(i)
            airy_argument = self._get_airy_argument(E, U, i)
            airy = special.airy(airy_argument)
            solution[0][i] = A * airy[0] + B * airy[2]
            solution[1][i] = A * airy[1] + B * airy[3]

        (A, B) = self._get_initial_right_AB(E, solution_start[2], solution_start[3])
        U = self._get_U(N - 1)
        for i in range(N - 1, N / 3, -1):
            if i < N - 1 and abs(self.potential[i] - self.potential[i+1]) > potential_threshold:
                (A, B) = self._get_AB(E, i + 1, i, solution[2][i + 1], solution[3][i + 1])
                U = self._get_U(i)
            airy_argument = self._get_airy_argument(E, U, i)
            airy = special.airy(airy_argument)
            solution[2][i] = A * airy[0] + B * airy[2]
            solution[3][i] = A * airy[1] + B * airy[3]
        return solution

    def _get_airy_argument(self, E, U, index):
        self._reset_to_default_units()
        x = self.length.value * index / (len(self.potential) - 1)
        E.convert_to(E.units_default)
        U.convert_to(E.units)
        return (x / self.x0.value) - (E.value - U.value) / self.eps0.value

    def _get_U(self, index):
        self._reset_to_default_units()
        potential_value = self.potential[index]
        x = self.length.value * index / (len(self.potential) - 1)
        U = spsc_data.EnergyValue(potential_value - self.slope.value * x * constants.e)
        return U

    def _get_initial_left_AB(self, E, func_value, der_value):
        return self._get_AB(E, 0, 0, func_value, der_value)

    def _get_initial_right_AB(self, E, func_value, der_value):
        index = len(self.potential) - 1
        U = self._get_U(index)
        airy_arg = self._get_airy_argument(E, U, index)
        airy = special.airy(airy_arg)
        A = func_value / airy[0]
        return A, 0

    def _get_AB(self, E, prev_index, index, func_value, der_value):
        U = self._get_U(index)
        airy_arg = self._get_airy_argument(E, U, prev_index)
        airy = special.airy(airy_arg)
        lin_matrix = np.array([[airy[0], airy[2]], [airy[1], airy[3]]])
        lin_right = np.array([[func_value], [der_value]])
        AB = np.linalg.solve(lin_matrix, lin_right)
        return AB[0], AB[1]


class SolutionIterationSlopePotentialFactory(ASolutionIterationFactory):
    def get_iteration(self, potential, mass, length):
        return SolutionIterationSlopePotential(potential, mass, length)
