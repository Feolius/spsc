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

# class MassUnitsConvertTestCase(unittest.TestCase):
#
#     def test_g_to_g(self):
#         units = spsc_units.MassUnits("g")
#         value = 2000
#         new_value = units.convert_value_to(value, "g")
#         self.assertEqual(new_value, 2000)
#
#     def test_kg_to_kg(self):
#         units = spsc_units.MassUnits("kg")
#         value = 2000
#         new_value = units.convert_value_to(value, "kg")
#         self.assertEqual(new_value, 2000)
#
#     def test_g_to_kg(self):
#         units = spsc_units.MassUnits("g")
#         value = 2
#         new_value = units.convert_value_to(value, "kg")
#         self.assertEqual(new_value, 2)
#
#     def test_kg_to_g(self):
#         units = spsc_units.MassUnits("kg")
#         value = 4
#         new_value = units.convert_value_to(value, "g")
#         self.assertEqual(new_value, 4000)