using MicroMagnetic
using CairoMakie
using DelimitedFiles

@using_gpu()

mesh =  FDMesh(nx=200, ny=50, nz=1, dx=2.5e-9, dy=2.5e-9, dz=3e-9);

sim = create_sim(mesh, name="std4", driver="SD", Ms=8.0e5, A=1.3e-11, demag=true,  m0=(1, 0.25, 0.1));
relax(sim, maxsteps=5000, stopping_dmdt=0.01)

plot_m(sim, component='x')

set_driver(sim, driver="LLG", alpha=0.02, gamma = 2.211e5)
add_zeeman(sim, (-24.6mT, 4.3mT, 0))
run_sim(sim, steps=100, dt=1e-11)

jld2movie("std4.jld2", output="assets/std4.mp4", component='x')

function plot_m_ts()
    folder = @__DIR__
    data = readdlm("std4_llg.txt", skipstart=2)
    oommf = readdlm("assets/std4_oommf.txt")

    fig = Figure(size = (800, 480))
    ax = Axis(fig[1, 1],
        xlabel = "Time (ns)",
        ylabel = "m"
    )

    lines!(ax, oommf[:,1]*1e9, oommf[:,2], label="OOMMF")
    lines!(ax, oommf[:,1]*1e9, oommf[:,3])
    lines!(ax, oommf[:,1]*1e9, oommf[:,4])

    scatter!(ax, data[:,2]*1e9, data[:,4], markersize = 6, label="MicroMagnetic")
    scatter!(ax, data[:,2]*1e9, data[:,5], markersize = 6)
    scatter!(ax, data[:,2]*1e9, data[:,6], markersize = 6)

    axislegend()

    return fig

end

plot_m_ts()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
