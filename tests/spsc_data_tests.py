import unittest
import src.spsc_data as spsc_data
import numpy as np


class PhysValueInstantiationTestCase(unittest.TestCase):
    class TestValue(spsc_data.APhysValue):

        @property
        def units_default(self):
            return "units_default"

        @property
        def units_options(self):
            return {"units_1": 100.0, "units_2": 0.1}

    def test_usual_instantiation(self):
        value = self.TestValue(5, "units_2")
        self.assertEqual(value.value, 5)
        self.assertEqual(value.units, "units_2")

    def test_default_units_instantiation(self):
        value = self.TestValue(3)
        self.assertEqual(value.value, 3)
        self.assertEqual(value.units_default, 'units_default')

    def test_wrong_units_instantiation(self):
        with self.assertRaises(ValueError) as context:
            value = self.TestValue(4, "units_3")


class PhysValueOperatorTestCase(unittest.TestCase):
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

    def test_add_same_units(self):
        value1 = self.TestValue(5)
        value2 = self.TestValue(4)
        value = value1 + value2
        self.assertEqual(value.value, 9)

    def test_add_different_units(self):
        value1 = self.TestValue(5)
        value2 = self.TestValue(4, "units_1")
        value = value1 + value2
        self.assertEqual(value.value, 45)

    def test_add_different_types(self):
        value1 = self.TestValue(5)
        value2 = self.AnotherTestValue(5)
        with self.assertRaises(ValueError) as context:
            value = value1 + value2

    def test_sub_same_units(self):
        value1 = self.TestValue(5)
        value2 = self.TestValue(4)
        value = value1 - value2
        self.assertEqual(value.value, 1)

    def test_sub_different_units(self):
        value1 = self.TestValue(4, "units_1")
        value2 = self.TestValue(5)
        value = value1 - value2
        self.assertEqual(value.value, 3.5)

    def test_sub_different_types(self):
        value1 = self.TestValue(5)
        value2 = self.AnotherTestValue(5)
        with self.assertRaises(ValueError) as context:
            value = value1 - value2

    def test_multiply_int(self):
        value = self.TestValue(5)
        value = value * 4
        self.assertEqual(value.value, 20)

    def test_multiply_float(self):
        value = self.TestValue(4)
        value = value * 5.5
        self.assertEqual(value.value, 22.0)

    def test_multiply_str(self):
        value = self.TestValue(4)
        with self.assertRaises(ValueError) as context:
            value = value * "str"

    def test_str(self):
        value = self.TestValue(4, "units_2")
        value_str = str(value)
        self.assertEqual(value_str, "4 units_2")


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

    class ValueArray(spsc_data.APhysValueArray, spsc_data.EmptyUnitsValue):
        pass

    def test_list_instantiation(self):
        arr = self.ValueArray([1, 2, 5])
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 2)
        self.assertEqual(arr.value[2], 5)

    def test_ndarray_instantiation(self):
        ndarr = np.array([1, 3, 5])
        arr = self.ValueArray(ndarr)
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 3)
        self.assertEqual(arr.value[2], 5)

    def test_ndarray_2_3_instantiation(self):
        ndarr = np.array([[1, 2, 3], [4, 5, 6]])
        arr = self.ValueArray(ndarr)
        self.assertEqual(arr.value[0], 1)
        self.assertEqual(arr.value[1], 2)
        self.assertEqual(arr.value[2], 3)
        self.assertEqual(arr.value[3], 4)
        self.assertEqual(arr.value[4], 5)
        self.assertEqual(arr.value[5], 6)



class PhysValueArrayOperatorTestCase(unittest.TestCase):

    class TestValue(spsc_data.APhysValue):

        @property
        def units_default(self):
            return "units_default"

        @property
        def units_options(self):
            return {"units_1": 10.0, "units_2": 0.1}

    class ValueArray(spsc_data.APhysValueArray, TestValue):
        pass

    class AnotherValueArray(spsc_data.APhysValueArray, TestValue):
        pass

    def test_get_item(self):
        arr = self.ValueArray([1, 2, 5])
        self.assertEqual(arr[0], 1)
        self.assertEqual(arr[1], 2)
        self.assertEqual(arr[2], 5)

    def test_in(self):
        arr = self.ValueArray([1, 2, 5])
        self.assertIn(1, arr)
        self.assertNotIn(3, arr)

    def test_len(self):
        arr = self.ValueArray([1, 2, 5])
        self.assertEqual(3, len(arr))

    def test_empty_len(self):
        arr = self.ValueArray([])
        self.assertEqual(0, len(arr))

    def test_add_same_units(self):
        arr1 = self.ValueArray([1, 2, 3])
        arr2 = self.ValueArray([3, 2, 1])
        arr = arr1 + arr2
        self.assertEqual(arr[0], 4)
        self.assertEqual(arr[1], 4)
        self.assertEqual(arr[2], 4)

    def test_add_different_units(self):
        arr1 = self.ValueArray([1, 2, 3])
        arr2 = self.ValueArray([3, 2, 1], "units_1")
        arr = arr1 + arr2
        self.assertEqual(arr[0], 31)
        self.assertEqual(arr[1], 22)
        self.assertEqual(arr[2], 13)

    def test_add_different_types(self):
        arr1 = self.ValueArray([1, 2, 3])
        arr2 = self.AnotherValueArray([3, 2, 1], "units_1")
        with self.assertRaises(ValueError) as context:
            arr = arr1 + arr2

    def test_sub_same_units(self):
        arr1 = self.ValueArray([1, 2, 3])
        arr2 = self.ValueArray([3, 2.5, 1])
        arr = arr1 - arr2
        self.assertEqual(arr[0], -2)
        self.assertEqual(arr[1], -0.5)
        self.assertEqual(arr[2], 2)

    def test_sub_different_units(self):
        arr1 = self.ValueArray([1, 2.5, 3])
        arr2 = self.ValueArray([3, 2, 1], "units_1")
        arr = arr1 - arr2
        self.assertEqual(arr[0], -29)
        self.assertEqual(arr[1], -17.5)
        self.assertEqual(arr[2], -7)

    def test_sub_different_types(self):
        arr1 = self.ValueArray([1, 2, 3])
        arr2 = self.AnotherValueArray([3, 2, 1], "units_1")
        with self.assertRaises(ValueError) as context:
            arr = arr1 - arr2

    def test_multiply_int(self):
        arr = self.ValueArray([1, -2, 3])
        arr = arr * 4
        self.assertEqual(arr[0], 4)
        self.assertEqual(arr[1], -8)
        self.assertEqual(arr[2], 12)

    def test_multiply_float(self):
        arr = self.ValueArray([1, -2, 3])
        arr = arr * 2.5
        self.assertEqual(arr[0], 2.5)
        self.assertEqual(arr[1], -5)
        self.assertEqual(arr[2], 7.5)

class PhysValueArrayConvertUnitsTestCase(unittest.TestCase):

    class TestValue(spsc_data.APhysValue):

        @property
        def units_default(self):
            return "units_default"

        @property
        def units_options(self):
            return {"units_1": 10.0, "units_2": 0.1}


    class TestValueArray(spsc_data.APhysValueArray, TestValue):
        pass

    def test_convert_1(self):
        arr = self.TestValueArray([1, 4, 5])
        arr.convert_to("units_2")
        self.assertEqual(arr[0], 10)
        self.assertEqual(arr[1], 40)
        self.assertEqual(arr[2], 50)

    def test_convert_2(self):
        arr = self.TestValueArray([1, 4, 5])
        arr.convert_to("units_1")
        self.assertEqual(arr[0], 0.1)
        self.assertEqual(arr[1], 0.4)
        self.assertEqual(arr[2], 0.5)
