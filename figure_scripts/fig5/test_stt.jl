using MicroMagnetic
using Printf
using DelimitedFiles
using CairoMakie

function init_dw(i,j,k,dx,dy,dz)
	if i < 250
		return (1,0.1,0)
  	elseif i< 260
   		return (0,1,0.1)
  	else
   		return (-1,0.1,0)
  	end
end

function run(;b=0.5e9, omega=2e9, alpha=0.1, beta=0.2, ux=100)

	function ufun(t)
		return exp(-b*t)*cos(omega*t)
	end

	mesh =  FDMesh(nx=500, ny=1, nz=1, dx=1e-9, dy=2e-9, dz=2e-9)
	sim = create_sim(mesh,  Ms=8e5, A=1.3e-11, Ku=1e5, axis=(1,0,0), m0=init_dw, name="dw")
	relax(sim, maxsteps=2000, stopping_dmdt=0.01, save_vtk_every = -1)

	set_driver(sim, driver="LLG_STT", alpha=alpha, beta=beta, ux=ux, ufun=ufun)

	function call_back_fun(sim, t)
    	m = Array(sim.spin)
		m = reshape(m, 3, 500)
		mx = m[1, :]
		id = argmin(abs.(mx))
		phi = atan(m[3, id],m[2, id])
    	open(@sprintf("phi_%g.txt", beta), "a") do f
        	write(f, @sprintf("%g  %g  %g\n", t, id, phi))
    	end
	end

	run_sim(sim, steps=50, dt=2e-10, save_m_every = -1, call_back=call_back_fun)
end


isfile("phi_0.txt") || run(beta=0.0)
isfile("phi_0.1.txt") || run(beta=0.1)
isfile("phi_0.2.txt") || run(beta=0.2)

function analytical_phi(beta, ts; b=0.5e9, w=1e9, alpha=0.1, ux=100)
	Delta = sqrt(1.3e-11/1e5)
    numerator = exp.(-b * ts) .* (-b * cos.(ts * w) + b * exp.(b * ts) + w * sin.(ts * w))
    denominator = (b^2 + w^2) * (1 + alpha^2) * Delta
    return numerator*(alpha-beta)*ux / denominator
end

function plot_phi()

    fig = Figure(size = (500, 360), fontsize = 18)
    ax = Axis(fig[1, 1],
        xlabel = "Time (ns)",
        ylabel = L"\phi"
    )

    ts = range(0, 10e-9, 501)
	
    phi0 = analytical_phi(0, ts, w=2e9)
	phi1 = analytical_phi(0.1, ts, w=2e9)
	phi2 = analytical_phi(0.2, ts, w=2e9)

    ts_ns = ts*1e9
    l1 = lines!(ax, linestyle = :solid, ts_ns, phi0, label=L"\beta=0")
    l2 = lines!(ax, linestyle = :dash,  linewidth=2, ts_ns, phi1, color=:sienna1, label=L"\beta=0.1")
    l3 = lines!(ax, linestyle = Linestyle([0.5, 1.0, 1.5, 2.5]), linewidth=2, color=:slateblue1, ts_ns, phi2, label=L"\beta=0.2")
    
	data0 = readdlm("phi_0.txt")
	data1 = readdlm("phi_0.1.txt")
	data2 = readdlm("phi_0.2.txt")
    s1 = scatter!(ax, data0[:, 1]*1e9, data0[:, 3] .- data0[1, 3], markersize = 8, label=L"\text{MicroMagnetic.jl}")
    scatter!(ax, data1[:, 1]*1e9, data1[:, 3] .- data1[1, 3], markersize = 8, color=:sienna1 )
    scatter!(ax, data2[:, 1]*1e9, data2[:, 3] .- data2[1, 3], markersize = 8, color=:slateblue1)

	axislegend(ax, [l1, l2, l3], [L"\beta=0", L"\beta=0.1", L"\beta=0.2"], position = :rt, orientation = :horizontal, labelsize=18)
	axislegend(ax, [s1], [L"\text{MicroMagnetic.jl}"], position=(0.75, 0.07), orientation = :horizontal, labelsize=16)

    save("dw_stt.pdf", fig)

end

plot_phi()
