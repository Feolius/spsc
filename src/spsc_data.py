from abc import ABCMeta, abstractproperty
import numpy as np
import matplotlib.pyplot as plt
import spsc_io


class APhysValue(spsc_io.Default):
    __metaclass__ = ABCMeta

    @abstractproperty
    def units_default(self):
        return ""

    @abstractproperty
    def units_options(self):
        return {}

    def __init__(self, value, units=None):
        if units is None:
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
                new_value = new_value * (1 / self.units_options[new_units])
            self.value = new_value
            self.units = new_units

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

    def __str__(self):
        return str(self.value) + " " + self.units

    def __eq__(self, other):
        equal = False
        if type(other) == self.__class__:
            equal = self.value == other.value and self.units == other.units
        return equal

    def __ne__(self, other):
        return not self == other

    def to_dict(self):
        dct = {
            "value": self.value,
            "units": self.units
        }
        return dct

    @classmethod
    def from_dict(cls, dct):
        if "value" not in dct:
            raise KeyError("Value key is not found in import dictionary")
        if "units" in dct:
            obj = cls(dct["value"], dct["units"])
        else:
            obj = cls(dct["value"])
        return obj

    @classmethod
    def get_all_units(cls):
        units = []
        obj = cls(0)
        units.append(obj.units_default)
        units.extend(obj.units_options.keys())
        return units


class EmptyUnitsValue(APhysValue):
    @property
    def units_default(self):
        return ""

    @property
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
        return {"eV": (1 / (6.24150965 * 10 ** 11)), "J": float(10 ** 7)}


class DensityValue(APhysValue):
    @property
    def units_default(self):
        return "cm^-2"

    @property
    def units_options(self):
        return {"m^-2": float(100 ** -2)}


class ChargeValue(APhysValue):
    @property
    def units_default(self):
        return "esu"

    @property
    def units_options(self):
        return {"C": 2997924580.0}


class APhysValueArray(APhysValue):
    __metaclass__ = ABCMeta

    def __init__(self, array, units=None):
        if type(array) == list:
            self.value = np.array(array, "float64")
        elif type(array) == np.ndarray:
            self.value = array.reshape(array.size).astype("float64")
        else:
            raise ValueError("Bad argument is passed when creating Value Array")
        super(APhysValueArray, self).__init__(self.value, units)

    def instant_plot(self):
        plt.gcf().clear()
        plt.plot(self.value)
        plt.ion()
        plt.show()

    def __getitem__(self, item):
        return self.value[item]

    def __setitem__(self, key, value):
        self.value[key] = value

    def __len__(self):
        return self.value.size

    def __contains__(self, item):
        return item in self.value

    def to_dict(self):
        dct = {
            "value": self.value.tolist(),
            "units": self.units
        }
        return dct

    def __eq__(self, other):
        equal = False
        if type(other) == self.__class__:
            equal = np.array_equal(self.value, other.value) and self.units == other.units
        return equal

    def __ne__(self, other):
        return not self == other

    def __delitem__(self, key):
        self.value = np.delete(self.value, key)

    def append(self, other):
        if type(other) == self.__class__:
            if other.units != self.units:
                other.convert_to(self.units)
            self.value = np.append(self.value, other.value)
        else:
            raise ValueError("Invalid argument type given for Value Array append method")

    @classmethod
    def from_dict(cls, dct):
        if "value" not in dct:
            raise KeyError("Value key is not found in import dictionary")
        if "units" in dct:
            obj = cls(dct["value"], dct["units"])
        else:
            obj = cls(dct["value"])
        return obj

    @classmethod
    def get_all_units(cls):
        units = []
        obj = cls([])
        units.append(obj.units_default)
        units.extend(obj.units_options.keys())
        return units


class Potential(APhysValueArray, EnergyValue):
    pass


class WaveFunction(APhysValueArray, EmptyUnitsValue):
    pass


class Density(APhysValueArray, DensityValue):
    pass

