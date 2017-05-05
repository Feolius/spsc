import spsc_data


class State(object):

    def __init__(self):
        self.electron_states = []
        self.static_density = spsc_data.Density([])
        self.density_potential = spsc_data.Potential([])
        self.length = 0

    @classmethod
    def from_dict(cls, dict):
        state = cls()
        if "electron_states" in dict:
            pass
        if "static_density" in dict:
            state.static_density = spsc_data.Density(dict["static_density"])
        if "density_potential" in dict:
            state.density_potential = spsc_data.Potential(dict["density_potential"])
        if "length" in dict:
            state.length = spsc_data.LengthValue()

