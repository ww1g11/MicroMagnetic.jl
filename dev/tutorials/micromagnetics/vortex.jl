using MicroMag
using CairoMakie

mesh =  FDMesh(dx=2e-9, dy=2e-9, dz=5e-9, nx=100, ny=100, nz=4);

geo = Cylinder(radius=100e-9);

function init_fun(i,j,k,dx,dy,dz)
  x = i-50.5
  y = j-50.5
  r = (x^2+y^2)^0.5
  if r<20
    return (0,0,1)
  end
  return (y/r, -x/r, 0)
end

sim = create_sim(mesh, shape=geo, Ms=8e5, A = 1.3e-11, demag=true, m0=init_fun);
relax(sim, maxsteps=5000, stopping_dmdt=0.1)
save_vtk(sim, "vortex")

plot_m(sim, figsize=(400, 400), arrows=(30,30))

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
