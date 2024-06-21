using MicroMagnetic
using Test
using NPZ

function analytical(alpha::Float64, gamma::Float64, H0::Float64, ts)
    precession = gamma / (1 + alpha * alpha)
    beta = precession * H0 * ts

    mx = cos.(beta) ./ cosh.(alpha .* beta)
    my = sin.(beta) ./ cosh.(alpha .* beta)
    mz = tanh.(alpha .* beta)
    return mx, my, mz
end

function test_llg_rk(;integrator="RungeKutta", steps=100, ts=1e-9)
    #Test mesh
    mesh = FDMesh(; nx=1, ny=1, dx=1e-9)

    sim = Sim(mesh; name="spin", integrator=integrator)
    set_Ms(sim, 8e5)
    sim.driver.alpha = 0.05
    sim.driver.gamma = 2.21e5
    sim.driver.integrator.dt = ts/steps

    add_zeeman(sim, (0, 0, 1e5))
    init_m0(sim, (1.0, 0, 0))

    for i in 1:steps
        MicroMagnetic.run_step(sim, sim.driver)
    end

    println(sim.driver.integrator.t == ts)

    #println(sim.spin[1]," ",sim.spin[2]," ",sim.spin[3])
    
    mx, my, mz = analytical(0.05, 2.21e5, 1e5, ts)
    m = Array(sim.spin)
    return ((m[1]-mx)^2 + (m[2]-my)^2 + (m[3]-mz)^2)^0.5
end

logrange(start,stepmul,length) = start .* stepmul .^ (0:(length-1))

if filesize("error.npy") == 0
    ts = 1e-9 
    Ns = logrange(20, 2, 10)
    errors = zeros(3, length(Ns))
    for (i, N) in enumerate(Ns)
        errors[1, i] = ts/N
        errors[2, i] = test_llg_rk(steps=N, integrator="RungeKutta")
        errors[3, i] = test_llg_rk(steps=N, integrator="RungeKuttaCayley")
        println(ts/N, " ", errors[:, i])
    end 
    npzwrite("error.npy", errors)

end

   

using CairoMakie

function plot_m()

    error = npzread("error.npy")

    fig = Figure(size = (500, 360), fontsize = 18)
    ax = Axis(fig[1, 1],
        xlabel = "Time step (s)",
        ylabel = "Absolute error",
        xscale = log10,
        yscale = log10
    )

    scatterlines!(ax, error[1,:], error[2,:], markersize = 12, marker = :rect, color = :sienna1, label="RungeKutta")
    scatterlines!(ax, error[1,:], error[3,:], markersize = 12, color=:slateblue1, strokecolor=:transparent, label="RungeKuttaCayley")

    axislegend(position=(0.95, 0.05), labelsize=16)

    save("error.pdf", fig)

end

plot_m()
