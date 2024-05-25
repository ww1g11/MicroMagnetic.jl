#stt dynamics script of fidimag as a comparison
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

import numpy as np
from fidimag.micro import Sim
from fidimag.common import CuboidMesh
from fidimag.micro import UniformExchange, Demag


mu0 = 4 * np.pi * 1e-7



def excite_system(mesh):

    # Specify the stt dynamics in the simulation
    sim = Sim(mesh, name='fidimag', driver='llg_stt')

    sim.driver.set_tols(rtol=1e-8, atol=1e-8)
    sim.driver.alpha = 0.1
    sim.driver.gamma = 2.211e5
    sim.Ms = 8.0e5
    sim.driver.p = 1

    sim.set_m(np.load('npys/m0_cpu.npy'))

    A = 1.3e-11
    exch = UniformExchange(A=A)
    sim.add(exch)

    sim.add(Demag())

    sim.driver.jx = -1e12
    sim.driver.beta = 0.05

    ts = np.linspace(0, 8e-9, 801)

    for t in ts:
        print('time', t)
        sim.driver.run_until(t)


if __name__ == '__main__':

    mesh = CuboidMesh(nx=20, ny=20, nz=2,dx=5, dy=5, dz=5.0, unit_length=1e-9)
    excite_system(mesh)
