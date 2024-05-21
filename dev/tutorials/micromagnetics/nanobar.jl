using MicroMag
using CairoMakie

mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=2.5e-9, nx=30, ny=5, nz=2);

sim = create_sim(mesh, driver="SD", name="bar", Ms=8e5, A=1e-12, demag=true, m0=(1, 1, 0));

plot_m(sim)

relax(sim, maxsteps=2000, stopping_dmdt=0.01)

fig = plot_m(sim)

save("bar.png", fig)

save_vtk(sim, "bar", fields=["exch", "demag"])

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
