import scipy.constants as constants
import spsc_data

h_plank = constants.hbar * (10 ** 7)
m_e = spsc_data.MassValue(constants.m_e, "kg")
m_e.convert_to("g")
m_e = m_e.value
e = spsc_data.ChargeValue(constants.elementary_charge, "C")
e.convert_to("esu")
e = e.value