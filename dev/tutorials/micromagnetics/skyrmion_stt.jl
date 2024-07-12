using MicroMagnetic
using Printf

@using_gpu()

mesh =  FDMesh(nx=400, ny=150, nz=1, dx=2e-9, dy=2e-9, dz=2e-9, pbc="xy");

sim = create_sim(mesh,  Ms=3.87e5, A=5.2e-12, D=1e-3, H=(0,0,160*mT), name="skx");

init_m0_skyrmion(sim, (200e-9, 80e-9), 2e-8)

relax(sim, maxsteps=20000, stopping_dmdt=0.01)

MicroMagnetic.save_vtk(sim, "skx")

set_driver(sim, driver="LLG_STT", alpha=0.05, beta=0.2, ux=-20)

function call_back_fun(sim, t)
    Rx, Ry = compute_guiding_center(sim)
    open("XY.txt", "a") do f
        write(f, @sprintf("%g  %g  %g\n", t, Rx, Ry))
    end
end

if !isfile("assets/skx.mp4")
    run_sim(sim, steps=100, dt=1e-10, save_m_every = 1, call_back=call_back_fun)
    jld2movie("skx.jld2", output="assets/skx.mp4")
end

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
