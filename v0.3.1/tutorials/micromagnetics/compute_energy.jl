using MicroMagnetic

mesh =  FDMesh(dx=5e-9, dy=5e-9, dz=5e-9, nx=10, ny=10, nz=1);

sim = create_sim(mesh, Ms=8.6e5, m0=(1,1,1));

demag = add_demag(sim);

zeeman = add_zeeman(sim, (0,0,1e5));

MicroMagnetic.effective_field(sim, sim.spin)

println("Demag Energy: ",sum(demag.energy), " J")
println("Zeeman Energy: ",sum(zeeman.energy), " J")

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl
