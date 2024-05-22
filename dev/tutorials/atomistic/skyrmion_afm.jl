using MicroMagnetic
using Printf

@using_gpu()

mesh =  CubicMesh(nx=800, ny=800, nz=1, dx=0.418e-9, dy=0.418e-9, dz=0.418e-9, pbc="xy")

function m0_fun(i, j, k, dx, dy, dz)
  flag = (i+j)%2 == 0 ? 1 : -1
  r2 = (i-300)^2 + (j-300)^2
  if r2 < 10^2
    return (0.1*flag, 0, -1*flag)
  end
  return (0.1,0,1*flag)
end

function relax_system()
  sim = Sim(mesh, driver="LLG")
  sim.driver.alpha = 0.5

  set_mu_s(sim, 2.96*mu_B)

  add_exch(sim, -34.4*meV)
  add_dmi(sim, 1.09*meV)
  add_anis(sim, 0.053*meV, axis=(0,0,1))
  #add_demag(sim)

  init_m0(sim, m0_fun)

  #using LLG to relax the system for the first 500 steps
  relax(sim, maxsteps=500, stopping_dmdt=1e-4)

  #change driver to SD since SD is much faster than LLG
  set_driver(sim, driver="SD")
  sim.driver.max_tau = 1
  relax(sim, maxsteps=50000, stopping_dmdt=1e-4)

  save_vtk(sim, "afm_skx")
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
