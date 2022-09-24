
"""
    set_mu_s(sim::AtomicSimGPU, Ms::NumberOrArrayOrFunction)

Set magnetic moment mu_s of the studied system. For example,

```julia
   set_mu_s(sim, 8.6e5)
```
or
```julia
function circular_shape(i,j,k,dx,dy,dz)
    if (i-50.5)^2 + (j-50.5)^2 <= 50^2
        return 8.6e5
    end
    return 0.0
end
set_mu_s(sim, circular_shape)
```
"""
function set_mu_s(sim::AtomicSimGPU, init::NumberOrArrayOrFunction)
    Float = _cuda_using_double.x ? Float64 : Float32
    Ms = zeros(Float, sim.nxyz)
    init_scalar!(Ms, sim.mesh, init)
    copyto!(sim.mu_s, Ms)
    return true
end

function set_mu_s_kagome(sim::AtomicSimGPU, Ms::Number)
    mesh = sim.mesh
    Float = _cuda_using_double.x ? Float64 : Float32
    mu_s = zeros(Float, sim.nxyz)
    for k = 1:mesh.nz, j = 1:mesh.ny, i = 1:mesh.nx
        id = index(i, j, k, mesh.nx, mesh.ny, mesh.nz)
        mu_s[id] = Ms
        if i%2==0 && j%2==0
            mu_s[id] = 0.0
        end
    end
    copyto!(sim.mu_s, mu_s)
    return true
end


"""
    add_exch(sim::AtomicSimGPU, J::Array; name="exch")

Add exchange energy to the system. The length of J should be equal to the length of neigbours.
"""
function add_exch(sim::AtomicSimGPU, J::Array; name="exch")
    nxyz = sim.nxyz
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*nxyz)
    energy = zeros(Float, nxyz)
    n_ngbs = sim.mesh.n_ngbs
    Js = CUDA.zeros(Float, n_ngbs)

    if length(J) != n_ngbs
        @error("The length of given Js is $(length(Js)) but we need an array with $a.")
    else
        copyto!(Js, [Float(i) for i in J])
    end

    exch = HeisenbergExchange(Js, field, energy, Float(0.0), name)
    push!(sim.interactions, exch)
    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return exch
end


"""
    add_exch(sim::AtomicSimGPU, J::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch(sim::AtomicSimGPU, J::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(n_ngbs)
    Js .= J
    add_exch(sim, Js, name=name)
end


"""
    add_dmi(sim::AtomicSimGPU, D::Real; name="dmi")

Add bulk dmi energy to the system.
"""
function add_dmi(sim::AtomicSimGPU, D::Real; name="dmi")
    nxyz = sim.nxyz
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*nxyz)
    energy = zeros(Float, nxyz)
    n_ngbs = sim.mesh.n_ngbs

    dmi = HeisenbergBulkDMI(Float(D), field, energy, Float(0.0), name)

    push!(sim.interactions, dmi)
    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
end


"""
    add_exch_kagome(sim::AtomicSimGPU, Jxy::Number, Jz::Number; name="exch")

Add exchange energy to the system.
"""
function add_exch_kagome(sim::AtomicSimGPU, Jxy::Number, Jz::Number; name="exch")
    n_ngbs = sim.mesh.n_ngbs
    Js = zeros(n_ngbs)

    if n_ngbs!=8
        error("The number of neigbours is not 8.")
    end

    Js[1:6] .= Jxy
    Js[7:8] .= Jz
    add_exch(sim, Js, name=name)
end

"""
    add_anis_kagome(sim::AtomicSimGPU, Ku::Float64; ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0), ax3=(-0.5,sqrt(3)/2,0), name="anis")

Add Anisotropy for kagome system, where the energy density is given by

```math
E_\\mathrm{anis} = - K_{u} (\\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis_kagome(sim::AtomicSimGPU, Ku::Float64; ax1=(-0.5,-sqrt(3)/2,0), ax2=(1,0,0), ax3=(-0.5,sqrt(3)/2,0), name="anis")
  nxyz = sim.nxyz
  T = _cuda_using_double.x ? Float64 : Float32
  field = zeros(T, 3*nxyz)
  energy = zeros(T, nxyz)
  anis =  KagomeAnisotropy(Float(Ku), ax1, ax2, ax3, field, energy, T(0.0), name)
  push!(sim.interactions, anis)

  if sim.save_data
      push!(sim.saver.headers, string("E_",name))
      push!(sim.saver.units, "J")
      id = length(sim.interactions)
      push!(sim.saver.results, o::AbstractSim->sum(o.interactions[id].energy))
  end
  return anis
end


function add_thermal_noise(sim::AtomicSimGPU, T::NumberOrArrayOrFunction; name="thermal", k_B=k_B)
    nxyz = sim.nxyz
    Float = _cuda_using_double.x ? Float64 : Float32
    field = zeros(Float, 3*nxyz)
    energy = zeros(Float, nxyz)
    Spatial_T = CUDA.zeros(Float, nxyz)
    eta = CUDA.zeros(Float, 3*nxyz)
    init_scalar!(Spatial_T , sim.mesh, T)
    thermal = StochasticFieldGPU(Spatial_T, eta, field, energy, Float(0.0), -1, name, k_B)
  
    push!(sim.interactions, thermal)
  
    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return thermal
end

"""
    add_magnetoelectric_laser(sim::AtomicSimGPU, lambda::Float64, E::Float64, B::Float64, omega::Float64; delta=0, direction=001, name="lasers")

Add the interaction of high-frequency lasers to multiferroic insulator Cu2OSeO3. 
The Hamiltonian is given by 

```math
\\mathcal{H}_\\mathrm{laser} =  -\\sum_{i} \\mu_s \\mathbf{m}_i \\cdot \\mathbf{B}(t) - \\sum_{i} \\mathbf{P}_i \\cdot \\mathbf{E}(t)
```

The high-frequency laser is described as 

```math
\\mathbf{E}(t) =  E ( \\sin (\\omega t + \\delta), \\cos \\omega t, 0) \\qquad
\\mathbf{B}(t) =  B ( \\cos \\omega t, -\\sin(\\omega t + \\delta), 0)
```
where δ determines the laser polarization, i.e., δ = 0 for right-circularly polarized (RCP), 
δ=π/2 for linearly polarized and δ=π for left-circularly polarized (LCP).
"""
function add_magnetoelectric_laser(sim::AtomicSimGPU, lambda::Float64, E::Float64, B::Float64, omega::Float64; delta=0, direction=001, name="lasers")
    nxyz = sim.nxyz
    F = _cuda_using_double.x ? Float64 : Float32
    field = zeros(F, 3*nxyz)
    energy = zeros(F, nxyz)

    laser = MagnetoelectricLaser(F(lambda), F(E),  F(B), F(omega), F(delta), 
                                direction, field, energy, F(0.0), name)
  
    push!(sim.interactions, laser)
  
    if sim.save_data
        push!(sim.saver.headers, string("E_",name))
        push!(sim.saver.units, "J")
        id = length(sim.interactions)
        push!(sim.saver.results, o::AbstractSim->o.interactions[id].total_energy)
    end
    return laser
end

