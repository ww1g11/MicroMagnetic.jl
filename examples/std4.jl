using SpinDynamics
using Printf
using NPZ

mesh =  FDMesh(nx=200, ny=50, nz=1, dx=2.5, dy=2.5, dz=3, unit_length=1e-9)

function relax_system()
  sim = Sim(mesh, name="std4_relax")
  set_Ms(sim, 8.0e5)
  sim.driver.alpha = 0.5
  sim.driver.precession = false

  add_exch(sim, 1.3e-11)
  add_demag(sim)

  init_m0(sim, (1, 0.25, 0.1))

  relax(sim, maxsteps=5000, stopping_dmdt = 0.1, save_m_every=1)
  npzwrite("m0.npy", sim.spin)
end

function run_dynamics()
  sim = Sim(mesh, name="std4_dyn")
  set_Ms(sim, 8.0e5)
  sim.driver.alpha = 0.02
  sim.driver.gamma = 2.211e5

  m0 = npzread("m0.npy")
  init_m0(sim, m0)
  add_exch(sim, 1.3e-11)
  add_demag(sim)

  mT = 0.001 / (4*pi*1e-7)
  add_zeeman(sim, (-24.6 * mT, 4.3 * mT, 0))

  for i = 1:100
    run_until(sim, 1e-11*i)
    println("step=",i)
  end
end


#relax_system()
run_dynamics()