using MicroMagnetic
using TimerOutputs
using Printf

@using_gpu()
MicroMagnetic.set_float(Float32)

function compute_time(;N=100, Nz=1)
    mesh =  FDMesh(nx=N, ny=N, nz=Nz, dx=4e-9, dy=4e-9, dz=4e-9);
    sim = create_sim(mesh, Ms=8e5, demag=true, A=1e-13, H=(0,100mT,0))
    set_driver(sim, driver="LLG", alpha=0.01, gamma = 2.211e5, integrator="Heun")

    #warm up
    for i=1:10
        MicroMagnetic.run_step(sim, sim.driver)
    end

    MicroMagnetic.reset_timer()

    # measure
    NN = N < 1000 ? 1000 : 500
    for i=1:NN
        @timeit MicroMagnetic.timer "run_step" MicroMagnetic.run_step(sim, sim.driver)
    end

    println(MicroMagnetic.timer)
    run_time = TimerOutputs.time(MicroMagnetic.timer["run_step"])/1e9 # seconds
    Ncalls = 2*TimerOutputs.ncalls(MicroMagnetic.timer["run_step"]) # two evals per step
    
    N2 = N*N
    return Ncalls, Ncalls*N2*Nz/run_time, run_time
end

f = open("cuda_time.txt","w")
write(f, "#N   Ncalls  throughput(cells/s)  total_time(s) \n")
for n in 6:-1:1
    N = 1000*n
    Ncalls, th, time = compute_time(N=N, Nz=1)
    GC.gc()
    write(f, @sprintf("%g  %g  %g  %g\n", N*N, Ncalls, th, time))
end
close(f)
