using CairoMakie
using DelimitedFiles
using CubicSplines
using MicroMagnetic

mesh = FDMesh(nx=60, ny=60, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy")
params = Dict(
    :Ms => 3.84e5,
    :A => 3.25e-12,
    :D => 5.83e-4,
    :H => (0, 0, 120 * mT)
)

function relax_skx()
    function m0_fun_skx(i, j, k, dx, dy, dz)
        r2 = (i - 30)^2 + (j - 30)^2
        if r2 < 10^2
            return (0.01, 0, -1)
        end
        return (0, 0, 1)
    end

    sim = create_sim(mesh; m0=m0_fun_skx, params...)
    relax(sim; maxsteps=2000, stopping_dmdt=0.01)
    save_vtk(sim, "skx")

    return plot_m(sim)
end

relax_skx()

init_images = [read_vtk("skx.vts"), (0, 0, 1)];

interpolation = [6];

sim = create_sim(mesh; params...);

neb = NEB(sim, init_images, interpolation; name="skx_fm", driver="SD");

relax(neb; stopping_dmdt=0.1, save_vtk_every=1000, maxsteps=5000)

function extract_data(; id=1)
    energy = readdlm("assets/skx_fm_energy.txt"; skipstart=2)
    dms = readdlm("assets/skx_fm_distance.txt"; skipstart=2)
    xs = zeros(length(dms[1, 1:end]))
    for i in 2:length(xs)
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

function plot_energy()
    fig = Figure(; resolution=(800, 480))
    ax = Axis(fig[1, 1]; xlabel="Distance (a.u.)", ylabel="Energy (meV)")

    xs, energy, xs2, energy2 = extract_data(; id=1)
    scatter!(ax, xs, energy; markersize=6, label="Initial energy path")
    lines!(ax, xs2, energy2)

    xs, energy, xs2, energy2 = extract_data(; id=500)
    scatter!(ax, xs, energy; markersize=6, label="Minimal energy path")
    lines!(ax, xs2, energy2)
    #linescatter!(ax, data[:,2]*1e9, data[:,5], markersize = 6)
    #linescatter!(ax, data[:,2]*1e9, data[:,6], markersize = 6)

    axislegend()

    save("energy.png", fig)

    return fig
end

plot_energy()

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
