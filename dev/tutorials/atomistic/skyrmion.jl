using MicroMag
using CairoMakie

mesh =  CubicMesh(nx=120, ny=120, nz=1, pbc="xy");

function m0_fun(i,j,k, dx, dy, dz)
  r2 = (i-60)^2 + (j-60)^2
  if r2 < 10^2
    return (0.1, 0, -1)
  end
  return (0,0,1)
end

function relax_system(mesh)
  #We create a simulation with 'SD' driver
  sim = Sim(mesh, driver="SD", name="skyrmion")

  #We set mu_s of the system
  set_mu_s(sim, 1.0)
  sim.driver.max_tau = 1000.0

  #Set the exchange, dmi and zeeman
  add_exch(sim, 1.0, name="exch")
  add_zeeman(sim, (0,0,3.75e-3))
  add_dmi(sim, 0.09, name="dmi")

  #Initialize the system using the `m0_fun` function
  init_m0(sim, m0_fun)

  #Relax the system
  relax(sim, maxsteps=2000, stopping_dmdt=1e-4, using_time_factor=false)

  return sim
end

sim = relax_system(mesh);

plot_m(sim)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
