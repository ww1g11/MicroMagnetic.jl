using MicroMagnetic
using Test

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg()

    mesh = FDMesh(nx=1, ny=1, nz=1)

    sim = Sim(mesh; name="spin")

    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.2
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.tol = 1e-8

    add_zeeman(sim, (0, 0, 1e5))

    init_m0(sim, (1.0, 0, 0))
    
    ts = range(0, 1e-9, 51)

    m = zeros(3, length(ts))
    for (i, t) in enumerate(ts)
        run_until(sim, t)
        m[:, i] .= sim.spin[:]
    end

    return ts, m
end


using CairoMakie

function plot_m()

    fig = Figure(size = (500, 360), fontsize = 18)
    ax = Axis(fig[1, 1],
        xlabel = "Time (ns)",
        ylabel = "m"
    )

    ts = range(0, 1e-9, 501)
    mx, my, mz = analytical(0.2, 2.21e5, 1e5, ts)

    ts_ns = ts*1e9
    l1 = lines!(ax, linestyle = :solid, ts_ns, mx, label=L"m_x")
    l2 = lines!(ax, linestyle = :dash,  linewidth=2, ts_ns, my, color=:sienna1, label=L"m_y")
    l3 = lines!(ax, linestyle = Linestyle([0.5, 1.0, 1.5, 2.5]), linewidth=2, color=:slateblue1, ts_ns, mz, label=L"m_z")
    
    ts, m = test_llg()
    ts_ns = ts*1e9
    s1 = scatter!(ax, ts_ns, m[1,:], markersize = 8, label=L"\text{MicroMagnetic.jl}")
    scatter!(ax, ts_ns, m[2,:], markersize = 8, color=:sienna1 )
    scatter!(ax, ts_ns, m[3,:], markersize = 8, color=:slateblue1)

    #axislegend(position=(0.95, 0.9), labelsize=18)
    axislegend(ax, [l1, l2, l3], [L"m_x", L"m_y", L"m_z"], position=(0.935, 0.87), orientation = :horizontal, labelsize=18)
	axislegend(ax, [s1], [L"\text{MicroMagnetic.jl}"], position=(0.84, 0.08), orientation = :horizontal, labelsize=16)

    save("m_ts.pdf", fig)

end

plot_m()
