import spsc_core
import spsc_data_generator

state = spsc_core.StateSimple.import_file("data/single_well26nm.yml")
solver = spsc_core.ShrodSolverSimple(state)
solver.solve()

