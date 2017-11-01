from abc import ABCMeta, abstractmethod
import spsc_data
import numpy as np
import spsc_constants as constants


class ASolutionStrategy(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def solve(self, density, eps, length):
        pass


class GaussSolutionStrategy(ASolutionStrategy):

    def solve(self, density, eps, length):
        N = len(density)
        h = 1.0 / (N - 1)
        length.convert_to(length.units_default)
        density.convert_to(density.units_default)
        potential = spsc_data.Potential(np.zeros((N,), "float64"))
        k = 2.0 * np.pi * (constants.e ** 2) * length.value * h / eps
        for i in xrange(N):
            for j in xrange(N):
                if i > j:
                    potential[i] += k * density[j] * (i - 2*j)
                else:
                    potential[i] += -k * density[j] * i
        return potential