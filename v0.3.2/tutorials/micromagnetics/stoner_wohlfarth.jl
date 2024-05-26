using MicroMagnetic

#We create a mesh for a cubic geometry 4nm x 4nm x 4nm
mesh = FDMesh(nx=4, ny=4, nz=4, dx=1e-9, dy=1e-9, dz=1e-9);

sim = create_sim(mesh, name="sw", driver="SD", Ms=1.0e6, A=1.3e-11, m0=(-1,1,0), Ku=5e4, axis=(1,1,0), H=(0,0,0))

#For each field, we relax the system to obtain its equilibrium state.
for i=-100:5:100
    Hx = i*mT # A/m
    update_zeeman(sim, (Hx,0,0))
    #Relax the system with stopping_dmdt=0.05, the write_data function will be called if save_m_every is positive
    relax(sim, maxsteps=10000, stopping_dmdt=0.05, save_data_every=-1)
end

using DelimitedFiles
using CairoMakie

function plot_loop()

    data = readdlm("./sw_sd.txt", skipstart=2)
    m, H = data[:, 3], data[:, 8]

    fig = Figure(size = (600, 400))
    ax = Axis(fig[1, 1],
        xlabel = "H (A/m)",
        ylabel = "mx"
    )

    scatterlines!(ax, H, m, markersize = 8, color = :blue, markercolor = :orange)
    scatterlines!(ax, -H, -m, markersize = 8, color = :blue, markercolor = :orange)

    expected = 39788.736 # A/m
    vlines!(ax, [expected, -expected], color = :red, linestyle = :dash)

    return fig

end

plot_loop()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
