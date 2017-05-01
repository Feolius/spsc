import spsc_data as spd

meters = spd.LengthUnits("m")
value = meters.convert_value_to(1, "cm")
print value