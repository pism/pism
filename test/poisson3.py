import PISM
import PISM.testing as pt
import numpy as np

config = PISM.Context().config

def run(Mx, Mz):
    L = 1.0
    grid = pt.shallow_grid(Mx=Mx, My=Mx, Lx=L, Ly=L)

    pp = PISM.Poisson3(grid, Mz, 1)

    inp = PISM.StressBalanceInputs()

    pp.update(inp, False)

    return pp

def write(model):

    output = PISM.util.prepare_output(config.get_string("output.file_name"))

    model.exact().write(output)
    model.solution().write(output)

    output.close()

if __name__ == "__main__":
    Mx = int(config.get_number("grid.Mx"))
    Mz = int(config.get_number("grid.Mz"))

    model = run(Mx, Mz)

    e = model.error()

    write(model)

    if model.grid().rank() == 0:
        print(f"Error = {e}")
