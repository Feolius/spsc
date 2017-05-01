from abc import ABCMeta, abstractmethod, abstractproperty
import numpy as np
import matplotlib.pyplot as plt


class APhysValue(object):
    __metaclass__ = ABCMeta

    @abstractproperty
    def units_default(self):
        return ""

    @abstractproperty
    def units_options(self):
        return {}

    def __init__(self, value, units=False):
        if not units:
            units = self.units_default
        if self._is_valid_units(units):
            self.units = units
            self.value = value
        else:
            raise ValueError("Unknown units value is passed during class instantiation.")

    def _is_valid_units(self, units):
        return (units == self.units_default) or (units in self.units_options.keys())

    def convert_to(self, new_units):
        if not self._is_valid_units(new_units):
            raise ValueError("Unknown units value is passed during units conversion.")
        new_value = self.value
        if new_units != self.units:
            units_default_coeff = 1
            if self.units != self.units_default:
                units_default_coeff = self.units_options[self.units]
            new_value *= units_default_coeff
            if new_units != self.units_default:
                new_value /= self.units_options[new_units]
            self.value = new_value

    def __add__(self, other):
        if type(other) == self.__class__:
            if other.units != self.units:
                other.convert_to(self.units)
            return self.__class__(self.value + other.value, self.units)
        else:
            raise ValueError("Invalid argument type given for + operator")

    def __sub__(self, other):
        if type(other) == self.__class__:
            if other.units != self.units:
                other.convert_to(self.units)
            return self.__class__(self.value - other.value, self.units)
        else:
            raise ValueError("Invalid argument type given for + operator")

    def __mul__(self, other):
        if type(other) == int or type(other) == float:
            return self.__class__(self.value * other, self.units)
        else:
            raise ValueError("Invalid argument type given for * operator")

class EmptyUnitsValue(APhysValue):

    def units_default(self):
        return ""

    def units_options(self):
        return {}


class LengthValue(APhysValue):
    @property
    def units_default(self):
        return "cm"

    @property
    def units_options(self):
        return {"m": 100.0, "mm": 0.1, "nm": float(10 ** -7)}


class MassValue(APhysValue):
    @property
    def units_default(self):
        return "g"

    @property
    def units_options(self):
        return {"kg": 1000.0}

class EnergyValue(APhysValue):
    @property
    def units_default(self):
        return "erg"

    @property
    def units_options(self):
        return {"eV": (1 / (6.24150965 * 10 ** 11)), "J" : float(10 ** 7)}

class APhysValueArray(APhysValue):

    def __init__(self, array, units=False):
        if type(array) == list:
            self.value = np.array(array).reshape((1, -1))
        elif type(array) == np.ndarray:
            self.value = array.reshape((1, -1))
        else:
            raise ValueError("Bad argument is passed when creating Value Array")
        super(APhysValueArray, self).__init__(array, units)

    def instant_plot(self):
        plt.plot(self.value)
        plt.show()

class Potential(APhysValueArray, EnergyValue):
    pass

