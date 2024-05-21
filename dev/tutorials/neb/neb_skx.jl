using Printf
using NPZ
using CairoMakie
using DelimitedFiles
using CubicSplines

using MicroMag
#MicroMag.cuda_using_double(true);

function create_sim(init_m_fun=(0,0,1))
    mesh =  FDMesh(nx=60, ny=60, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")
    sim = Sim(mesh, name="neb", driver="SD")
    set_Ms(sim, 3.84e5)

    init_m0(sim, init_m_fun)

    add_exch(sim, 3.25e-12)
    add_dmi(sim, 5.83e-4)
    add_zeeman(sim, (0, 0, 120*mT))

    return sim
end

function relax_skx()
    function m0_fun_skx(i,j,k, dx, dy, dz)
        r2 = (i-30)^2 + (j-30)^2
        if r2 < 10^2
          return (0.01, 0, -1)
        end
        return (0,0,1)
    end

    sim = create_sim(m0_fun_skx)
    relax(sim, maxsteps=2000, stopping_dmdt=0.01)
    npzwrite("skx.npy", Array(sim.spin))

    save_vtk(sim, "skx")
end

function plot_spatial_m(m; nx=60, ny=60, filename="")

  points = [Point3f(i, j, 0) for i in 1:2:nx for j in 1:2:ny]

  m = reshape(m, 3, nx, ny)
  mf = [Vec3f(m[1, i, j], m[2, i,j], m[3, i,j]) for i in 1:2:nx for j in 1:2:ny]
  mz = [m[3, i, j]  for i in 1:2:nx for j in 1:2:ny]

  fig = Figure(resolution = (800, 800))
  ax = Axis(fig[1, 1], backgroundcolor = "white")

  arrows!(ax, points, mf, fxaa=true, # turn on anti-aliasing
          color = vec(mz), linewidth = 0.5, arrowsize = 1, lengthscale = 1,
          align = :center
      )

  if length(filename)>0
    save(filename*".png", fig)
  end

  return fig

end

relax_skx()

plot_spatial_m(npzread("skx.npy"))

sim = create_sim()

init_images = [npzread("skx.npy"),  (0, 0, 1)]

interpolation  = [6]

if !isfile("skx_fm_energy.txt")
  #relax(neb, stopping_dmdt=0.1, save_vtk_every=1000, maxsteps=5000)
end

function extract_data(;id=1)
  energy = readdlm("assets/skx_fm_energy.txt", skipstart=2)
  dms = readdlm("assets/skx_fm_distance.txt", skipstart=2)
  xs = zeros(length(dms[1, 1:end]))
  for i=2:length(xs)
    xs[i] = sum(dms[id, 2:i])
  end

  et = energy[id, 2:end]
  e0 = minimum(et)
  energy_eV = (et .- e0) / meV

  spline = CubicSpline(xs, energy_eV)

  xs2 = range(xs[1], xs[end], 100)
  energy2 = spline[xs2]

  return xs, energy_eV, xs2, energy2
end


function plot_m()

  fig = Figure(resolution = (800, 480))
  ax = Axis(fig[1, 1],
      xlabel = "Distance (a.u.)",
      ylabel = "Energy (meV)"
  )

  xs, energy, xs2, energy2 = extract_data(id=1)
  scatter!(ax, xs, energy, markersize = 6, label="Initial energy path")
  lines!(ax, xs2, energy2)

  xs, energy, xs2, energy2 = extract_data(id=500)
  scatter!(ax, xs, energy, markersize = 6, label="Minimal energy path")
  lines!(ax, xs2, energy2)
  #linescatter!(ax, data[:,2]*1e9, data[:,5], markersize = 6)
  #linescatter!(ax, data[:,2]*1e9, data[:,6], markersize = 6)

  axislegend()

  save("energy.png", fig)

  return fig

end

plot_m()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
