using MicroMagnetic
using Printf
using SpecialFunctions
using CairoMakie
using StatsBase
using LinearAlgebra

@using_gpu()

function run()
	mesh =  FDMesh(nx=60, ny=60, nz=60, dx=2.5e-9, dy=3.2e-9, dz=3.5e-9)
	sim = create_sim(mesh, Ms=1.42e5, Ku=7.2e5, axis=(0,0,1), driver="LLG", alpha=0.1, T=300, name="sllg")
	run_until(sim, 1e-10)
    save_vtk(sim, "1.vts")
end

#run()

function compute_Z()
    K = 7.2e5
    V = 2.8e-26
    T = 300

    k_B = 1.3806488e-23
    chi =  K*V/(k_B*T)
    Z = 2*dawson(sqrt(chi))/sqrt(chi) 
    println(chi, "  ", Z)
    return chi, Z
end


function analytical()
    chi, Z = compute_Z()
    mzs = range(-1, 1, 201)
    ps = 1.0/Z*exp.(-chi *(1 .- mzs.^2))
    return mzs, ps
end

function plot_distribution()
    
    m = MicroMagnetic.read_vtk("1.vts")
    m = reshape(m, 3, div(length(m),3))

    hist = fit(Histogram, m[3, :], -1:0.1:1, closed=:right)
    mz = midpoints(hist.edges[1])
    h = normalize(hist, mode=:pdf)

    #with_theme(theme_latexfonts()) do
        fig = Figure(size = (500, 360), fontsize = 18)
        ax = Axis(fig[1, 1],
            xlabel = L"$m_z$",
            ylabel = L"log($P_\mathrm{eq}$)",
        )

        mzs, ps = analytical()
        l1 = lines!(ax, linestyle = :solid, mzs, log.(ps), color=:slateblue1, label="Analytical")
        s1 = scatter!(ax, mz, log.(h.weights), markersize = 10, strokewidth=1, alpha=0, color=:white, label="MicroMagnetic.jl")

        axislegend(ax, position = (0.5, 0.75), labelsize=14)
	    #axislegend(ax, [l1, l2, l3], [L"\beta=0", L"\beta=0.1", L"\beta=0.2"], position = :rt, orientation = :horizontal, labelsize=18)
	    #axislegend(ax, [s1], [L"\text{JuMag}"], position=(0.78, 0.07), orientation = :horizontal, labelsize=16)

        save("P_mz.pdf", fig)
    #end

end

plot_distribution()