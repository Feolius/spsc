import unittest
import src.spsc_data as spsc_data


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
        value1 = self.TestValue(5)
        value2 = self.TestValue(40, "units_1")
        value = value2 - value1
        self.assertEqual(value.value, 35)

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