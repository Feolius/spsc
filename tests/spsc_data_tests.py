import unittest
import src.spsc_data as spsc_data
import numpy as np
import random
import os


def random_list(length=None):
    if length is None:
        length = random.randrange(1, 10)
    rnd_lst = []
    for i in range(length):
        rnd_lst.append(random.uniform(0, 100))
    return rnd_lst


class TestValue(spsc_data.APhysValue):
    @property
    def units_default(self):
        return "units_default"

    @property
    def units_options(self):
        return {"units_1": 10.0, "units_2": 0.1}


class AnotherTestValue(spsc_data.APhysValue):
    @property
    def units_default(self):
        return "units_default"

    @property
    def units_options(self):
        return {"units_1": 100.0, "units_2": 0.1}


class TestValueArray(spsc_data.APhysValueArray, TestValue):
    pass


class AnotherTestValueArray(spsc_data.APhysValueArray, TestValue):
    pass


class PhysValueInstantiationTestCase(unittest.TestCase):

    def test_usual_instantiation(self):
        value = TestValue(5, "units_2")
        self.assertEqual(value.value, 5)
        self.assertEqual(value.units, "units_2")

    def test_default_units_instantiation(self):
        value = TestValue(3)
        self.assertEqual(value.value, 3)
        self.assertEqual(value.units_default, 'units_default')

    def test_wrong_units_instantiation(self):
        with self.assertRaises(ValueError) as context:
            value = TestValue(4, "units_3")

    def test_get_all_units(self):
        units = TestValue.get_all_units()
        test_units = ["units_default", "units_1", "units_2"]
        self.assertEqual(len(units), len(test_units))
        for i in range(len(units)):
            self.assertIn(units[i], test_units)


class PhysValueOperatorTestCase(unittest.TestCase):

    def test_add_same_units(self):
        value1 = TestValue(5)
        value2 = TestValue(4)
        value = value1 + value2
        self.assertEqual(value.value, 9)

    def test_add_different_units(self):
        value1 = TestValue(5)
        value2 = TestValue(4, "units_1")
        value = value1 + value2
        self.assertEqual(value.value, 45)

    def test_add_different_types(self):
        value1 = TestValue(5)
        value2 = AnotherTestValue(5)
        with self.assertRaises(ValueError) as context:
            value = value1 + value2

    def test_sub_same_units(self):
        value1 = TestValue(5)
        value2 = TestValue(4)
        value = value1 - value2
        self.assertEqual(value.value, 1)

    def test_sub_different_units(self):
        value1 = TestValue(4, "units_1")
        value2 = TestValue(5)
        value = value1 - value2
        self.assertEqual(value.value, 3.5)

    def test_sub_different_types(self):
        value1 = TestValue(5)
        value2 = AnotherTestValue(5)
        with self.assertRaises(ValueError) as context:
            value = value1 - value2

    def test_multiply_int(self):
        value = TestValue(5)
        value = value * 4
        self.assertEqual(value.value, 20)

    def test_multiply_float(self):
        value = TestValue(4)
        value = value * 5.5
        self.assertEqual(value.value, 22.0)

    def test_multiply_str(self):
        value = TestValue(4)
        with self.assertRaises(ValueError) as context:
            value = value * "str"

    def test_str(self):
        value = TestValue(4, "units_2")
        value_str = str(value)
        self.assertEqual(value_str, "4 units_2")

    def test_equal(self):
        value1 = TestValue(5, "units_1")
        value2 = TestValue(5, "units_1")
        value3 = TestValue(4, "units_1")
        value4 = TestValue(5, "units_2")
        value5 = TestValue(4, "units_2")
        self.assertEqual(value1, value2)
        self.assertNotEqual(value1, value3)
        self.assertNotEqual(value1, value4)
        self.assertNotEqual(value1, value5)


class LengthValueConvertTestCase(unittest.TestCase):

    def test_m_to_cm(self):
        value = spsc_data.LengthValue(5, "m")
        value.convert_to("cm")
        self.assertEqual(value.value, 500)

    def test_mm_to_cm(self):
        value = spsc_data.LengthValue(10, "mm")
        value.convert_to("cm")
        self.assertEqual(value.value, 1)

    def test_nm_to_cm(self):
        value = spsc_data.LengthValue(5, "nm")
        value.convert_to("cm")
        self.assertEqual(value.value, 5 * 10 ** -7)

    def test_m_to_nm(self):
        value = spsc_data.LengthValue(5, "m")
        value.convert_to("nm")
        self.assertEqual(value.value, 5 * 10 ** 9)

    def test_nm_to_m(self):
        value = spsc_data.LengthValue(5, "nm")
        value.convert_to("m")
        self.assertEqual(value.value, 5 * 10 ** -9)

    def test_cm_to_cm(self):
        value = spsc_data.LengthValue(3, "cm")
        value.convert_to("cm")
        self.assertEqual(value.value, 3)

    def test_m_to_m(self):
        value = spsc_data.LengthValue(5, "m")
        value.convert_to("m")
        self.assertEqual(value.value, 5)


class MassUnitsConvertTestCase(unittest.TestCase):

    def test_g_to_g(self):
        value = spsc_data.MassValue(5, "g")
        value.convert_to("g")

    def test_g_to_kg(self):
        value = spsc_data.MassValue(500, "g")
        value.convert_to("kg")
        self.assertEqual(value.value, 0.5)

    def test_kg_to_g(self):
        value = spsc_data.MassValue(5, "kg")
        value.convert_to("g")
        self.assertEqual(value.value, 5000)


class EnergyUnitsConvertTestCase(unittest.TestCase):

    def test_erg_to_ert(self):
        value = spsc_data.EnergyValue(5, "erg")
        value.convert_to("erg")
        self.assertEqual(value.value, 5)

    def test_erg_to_ev(self):
        value = spsc_data.EnergyValue(7, "erg")
        value.convert_to("eV")
        assert_value = value.value / 10 ** 11
        self.assertAlmostEqual(assert_value, 7 * 6.24150965)

    def test_j_to_erg(self):
        value = spsc_data.EnergyValue(8, "J")
        value.convert_to("erg")
        self.assertEqual(value.value, 8 * 10 ** 7)


class DensityUnitsConvertTestCase(unittest.TestCase):

    def test_cm_to_cm(self):
        value = spsc_data.DensityValue(5)
        value.convert_to("cm^-2")
        self.assertEqual(value.value, 5)

    def test_cm_to_m(self):
        value = spsc_data.DensityValue(5)
        value.convert_to("m^-2")
        self.assertEqual(value.value, 50000)

class ChargeUnitsConvertTestCase(unittest.TestCase):

    def test_esu_to_esu(self):
        value = spsc_data.ChargeValue(5)
        value.convert_to("esu")
        self.assertEqual(value.value, 5)

    def test_C_to_esu(self):
        value = spsc_data.ChargeValue(5, "C")
        value.convert_to("esu")
        self.assertEqual(value.value, 5 * 2997924580)


class PhysValueArrayInstantiationTestCase(unittest.TestCase):

    def test_list_instantiation(self):
        arr = TestValueArray([1, 2, 5])
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 2)
        self.assertEqual(arr.value[2], 5)

    def test_ndarray_instantiation(self):
        ndarr = np.array([1, 3, 5])
        arr = TestValueArray(ndarr)
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 3)
        self.assertEqual(arr.value[2], 5)

    def test_ndarray_2_3_instantiation(self):
        ndarr = np.array([[1, 2, 3], [4, 5, 6]])
        arr = TestValueArray(ndarr)
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 2)
        self.assertEqual(arr.value[2], 3)
        self.assertEqual(arr.value[3], 4)
        self.assertEqual(arr.value[4], 5)
        self.assertEqual(arr.value[5], 6)

    def test_get_all_units(self):
        units = TestValueArray.get_all_units()
        test_units = ["units_default", "units_1", "units_2"]
        self.assertEqual(len(units), len(test_units))
        for i in range(len(units)):
            self.assertIn(units[i], test_units)


class PhysValueArrayOperatorTestCase(unittest.TestCase):

    def test_get_item(self):
        arr = TestValueArray([1, 2, 5])
        self.assertEqual(arr[0], 1)
        self.assertEqual(arr[1], 2)
        self.assertEqual(arr[2], 5)

    def test_in(self):
        arr = TestValueArray([1, 2, 5])
        self.assertIn(1, arr)
        self.assertNotIn(3, arr)

    def test_len(self):
        arr = TestValueArray([1, 2, 5])
        self.assertEqual(3, len(arr))

    def test_empty_len(self):
        arr = TestValueArray([])
        self.assertEqual(0, len(arr))

    def test_add_same_units(self):
        arr1 = TestValueArray([1, 2, 3])
        arr2 = TestValueArray([3, 2, 1])
        arr = arr1 + arr2
        self.assertEqual(arr[0], 4)
        self.assertEqual(arr[1], 4)
        self.assertEqual(arr[2], 4)

    def test_add_different_units(self):
        arr1 = TestValueArray([1, 2, 3])
        arr2 = TestValueArray([3, 2, 1], "units_1")
        arr = arr1 + arr2
        self.assertEqual(arr[0], 31)
        self.assertEqual(arr[1], 22)
        self.assertEqual(arr[2], 13)

    def test_add_different_types(self):
        arr1 = TestValueArray([1, 2, 3])
        arr2 = AnotherTestValueArray([3, 2, 1], "units_1")
        with self.assertRaises(ValueError) as context:
            arr = arr1 + arr2

    def test_sub_same_units(self):
        arr1 = TestValueArray([1, 2, 3])
        arr2 = TestValueArray([3, 2.5, 1])
        arr = arr1 - arr2
        self.assertEqual(arr[0], -2)
        self.assertEqual(arr[1], -0.5)
        self.assertEqual(arr[2], 2)

    def test_sub_different_units(self):
        arr1 = TestValueArray([1, 2.5, 3])
        arr2 = TestValueArray([3, 2, 1], "units_1")
        arr = arr1 - arr2
        self.assertEqual(arr[0], -29)
        self.assertEqual(arr[1], -17.5)
        self.assertEqual(arr[2], -7)

    def test_sub_different_types(self):
        arr1 = TestValueArray([1, 2, 3])
        arr2 = AnotherTestValueArray([3, 2, 1], "units_1")
        with self.assertRaises(ValueError) as context:
            arr = arr1 - arr2

    def test_multiply_int(self):
        arr = TestValueArray([1, -2, 3])
        arr = arr * 4
        self.assertEqual(arr[0], 4)
        self.assertEqual(arr[1], -8)
        self.assertEqual(arr[2], 12)

    def test_multiply_float(self):
        arr = TestValueArray([1, -2, 3])
        arr = arr * 2.5
        self.assertEqual(arr[0], 2.5)
        self.assertEqual(arr[1], -5)
        self.assertEqual(arr[2], 7.5)

    def test_equal(self):
        arr1 = TestValueArray([1, 3, 4])
        arr2 = TestValueArray([1, 3, 4])
        arr3 = TestValueArray([1, 3, 4], "units_2")
        arr4 = TestValueArray([1, 3, 5])
        arr5 = TestValueArray([1, 3])
        self.assertEqual(arr1, arr2)
        self.assertNotEqual(arr1, arr3)
        self.assertNotEqual(arr1, arr4)
        self.assertNotEqual(arr1, arr5)

    def test_non_equal(self):
        arr1 = TestValueArray([1, 3, 4])
        arr2 = TestValueArray([1, 3])
        arr3 = TestValueArray([1, 2, 5])
        self.assertNotEqual(arr1, arr2)
        self.assertNotEqual(arr1, arr3)

    def test_delete(self):
        original_list = random_list()
        arr = TestValueArray(original_list)
        rnd_int = random.randrange(0, len(original_list))
        del arr[rnd_int]
        del original_list[rnd_int]
        new_arr = TestValueArray(original_list)
        self.assertEqual(arr, new_arr)


class PhysValueArrayConvertUnitsTestCase(unittest.TestCase):

    def test_convert_1(self):
        arr = TestValueArray([1, 4, 5])
        arr.convert_to("units_2")
        self.assertEqual(arr[0], 10)
        self.assertEqual(arr[1], 40)
        self.assertEqual(arr[2], 50)

    def test_convert_2(self):
        arr = TestValueArray([1, 4, 5])
        arr.convert_to("units_1")
        self.assertEqual(arr[0], 0.1)
        self.assertEqual(arr[1], 0.4)
        self.assertEqual(arr[2], 0.5)


class PhysValueIOTestCase(unittest.TestCase):

    def test_to_dict_int(self):
        rnd = random.randint(0, 100)
        value = TestValue(rnd)
        dct = value.to_dict()
        self.assertEqual(dct["value"], rnd)
        self.assertEqual(dct["units"], "units_default")

    def test_to_dict_float(self):
        rnd = random.uniform(0, 100)
        value = TestValue(rnd, "units_1")
        dct = value.to_dict()
        self.assertEqual(dct["value"], rnd)
        self.assertEqual(dct["units"], "units_1")

    def test_from_dict_int(self):
        rnd = random.randint(0, 100)
        dct = {
            "value": rnd,
            "units": "units_1"
        }
        value = TestValue.from_dict(dct)
        self.assertEqual(value.value, rnd)
        self.assertEqual(value.units, "units_1")

    def test_from_dict_float(self):
        rnd = random.uniform(0, 100)
        dct = {
            "value": rnd,
            "units": "units_1"
        }
        value = TestValue.from_dict(dct)
        self.assertEqual(value.value, rnd)
        self.assertEqual(value.units, "units_1")

    def test_import_export(self):
        value = TestValue(random.uniform(0, 100))
        test_file_name = "test.yml"
        value.export_file(test_file_name)
        new_value = TestValue.import_file(test_file_name)
        self.assertEqual(value, new_value)
        os.remove(test_file_name)


class PhysValueArrayIOTestCase(unittest.TestCase):

    def test_to_dict(self):
        rnd_lst = random_list()
        arr = TestValueArray(rnd_lst)
        dct = arr.to_dict()
        self.assertEqual(len(arr), len(dct["value"]))
        self.assertEqual(arr.units, dct["units"])
        for i in range(len(dct["value"])):
            self.assertEqual(arr[i], dct["value"][i])

    def test_from_dict(self):
        rnd_lst = random_list()
        dct = {
            "value": rnd_lst,
            "units": "units_2"
        }
        arr = TestValueArray.from_dict(dct)
        self.assertEqual(arr.units, dct["units"])
        self.assertEqual(len(arr), len(dct["value"]))
        for i in range(len(dct["value"])):
            self.assertEqual(arr[i], dct["value"][i])

    def test_import_export(self):
        rnd_lst = random_list()
        arr = TestValueArray(rnd_lst)
        test_file = "test2.yml"
        arr.export_file(test_file)
        new_arr = TestValueArray.import_file(test_file)
        self.assertEqual(arr, new_arr)
        os.remove(test_file)
