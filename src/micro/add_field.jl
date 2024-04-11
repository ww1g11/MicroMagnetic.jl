export add_zeeman, update_zeeman, add_anis, update_anis, add_cubic_anis, add_exch, add_dmi, add_demag

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Add a static Zeeman energy to the simulation.
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    n_total = sim.n_total

    T = single_precision.x ? Float32 : Float64
    field = zeros(T, 3 * n_total)

    init_vector!(field, sim.mesh, H0)

    field_kb = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy_kb = KernelAbstractions.zeros(backend[], T, n_total)

    copyto!(field_kb, field)

    zeeman = Zeeman(field_kb, energy_kb, name)
    push!(sim.interactions, zeeman)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple) && length(H0) == 3
            field_item = SaverItem((string(name, "_Hx"), string(name, "_Hy"),
                                    string(name, "_Hz")), ("A/m", "A/m", "A/m"),
                                   o::AbstractSim -> H0)
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end

    @info "Static Zeeman has been added."

    return zeeman
end

"""
    update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")

Set the Zeeman field to H0 where H0 is TupleOrArrayOrFunction according to its name. For example,

```julia
   add_zeeman(sim, (0,0,0), name="my_H")  #create a zeeman energy with field (0,0,0) A/m
   update_zeeman(sim, (0,0,1e5), name="my_H")  #change the field to (0,0,1e5) A/m
```
"""
function update_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction; name="zeeman")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64
    field = zeros(T, 3 * n_total)
    init_vector!(field, sim.mesh, H0)

    for i in sim.interactions
        if i.name == name
            copyto!(i.field, field)
            return nothing
        end
    end
    return nothing
end

"""
    add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function; name="timezeeman")

Add a time varying zeeman to system.

The input `ft` is a function of time `t` and its return value should be a tuple with length 3.

Example:

```julia
  function time_fun(t)
    w = 2*pi*2.0e9
    return (sin(w*t), cos(w*t), 0)
  end

  function spatial_H(i, j, k, dx, dy, dz)
    H = 1e3
    if i<=2
        return (H, H, 0)
    end
    return (0, 0, 0)
  end

  add_zeeman(sim, spatial_H, time_fun)
```
"""
function add_zeeman(sim::AbstractSim, H0::TupleOrArrayOrFunction, ft::Function;
                    name="timezeeman")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64

    init_field = zeros(T, 3 * n_total)
    init_vector!(init_field, sim.mesh, H0)

    field_kb = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    field = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(backend[], T, n_total)

    copyto!(field_kb, init_field)

    zeeman = TimeZeeman(ft, field_kb, field, energy, name)
    push!(sim.interactions, zeeman)

    if sim.save_data
        id = length(sim.interactions)
        if isa(H0, Tuple) && length(H0) == 3  # FIXME: the output should depends on time!!!
            field_item = SaverItem((string(name, "_Hx"), string(name, "_Hy"),
                                    string(name, "_Hz")), ("A/m", "A/m", "A/m"),
                                   o::AbstractSim -> H0)
            push!(sim.saver.items, field_item)
        end
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return zeeman
end

"""
    add_exch(sim::AbstractSim, A::NumberOrTupleOrArrayOrFunction; name="exch")

Add exchange energy to the system.

# Examples:
```julia
    add_exch(sim, 1e-11)
```

or 

```julia
    add_exch(sim, (2e-12,5e-12,0))
```
"""
function add_exch(sim::AbstractSim, A::NumberOrTupleOrArrayOrFunction; name="exch")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64
    field = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(backend[], T, n_total)

    exch = nothing
    if isa(A, Number)
        exch = VectorExchange(T(A), T(A), T(A), field, energy, name)
    elseif isa(A, Tuple) && length(A) == 3
        exch = VectorExchange(T(A[1]), T(A[2]), T(A[3]), field, energy, name)
    else
        Spatial_A = zeros(T, sim.n_total)
        init_scalar!(Spatial_A, sim.mesh, A)

        A_kb = KernelAbstractions.zeros(backend[], T, n_total)
        copyto!(A_kb, Spatial_A)

        exch = Exchange(A_kb, field, energy, name)
    end

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    @info "Exchange has been added."
    return exch
end

@doc raw"""
    add_dmi(sim::AbstractSim, D::NumberOrTupleOrArrayOrFunction; name="dmi", type="bulk")

Add DMI to the system. `type` could be "bulk" or "interfacial"

The energy of interlayer DMI is defined as 

```math
E_\mathrm{dmi-int} =  \int_\Gamma \mathbf{D} \cdot \left(\mathbf{m}_{i} \times \mathbf{m}_{j}\right) dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$. 
$\mathbf{D}$ is the effective DMI vector. 

The effective field is given
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s \Delta}  \mathbf{D} \times \mathbf{m}_{j} 
```

Examples:

```julia
   add_dmi(sim, 1e-3, type="interfacial")
```
or
```julia
   add_dmi(sim, 1e-3, type="bulk")
```

```julia
   add_dmi(sim, (1e-3, 1e-3, 0), type="bulk")
```

"""
function add_dmi(sim::AbstractSim, D::NumberOrTupleOrArrayOrFunction; name="dmi",
                 type="bulk")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64
    field = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(backend[], T, n_total)

    if type == "bulk"
        if isa(D, Number)
            dmi = BulkDMI(T(D), T(D), T(D), field, energy, name)
        elseif isa(D, Tuple) && length(D) == 3
            dmi = BulkDMI(T(D[1]), T(D[2]), T(D[3]), field, energy, name)
        else
            Spatial_D = zeros(T, sim.n_total)
            init_scalar!(Spatial_D, sim.mesh, D)

            D_kb = KernelAbstractions.zeros(backend[], T, n_total)
            copyto!(D_kb, Spatial_D)

            dmi = SpatialBulkDMI(D_kb, field, energy, name)
        end

        @info "Bulk DMI has been added."
    elseif type == "interfacial"
        Spatial_D = zeros(T, sim.n_total)
        init_scalar!(Spatial_D, sim.mesh, D)

        D_kb = KernelAbstractions.zeros(backend[], T, n_total)
        copyto!(D_kb, Spatial_D)

        dmi = InterfacialDMI(D_kb, field, energy, name)
        @info "Interfacial DMI has been added."

    else
        error("Supported DMI type:", "interfacial", "bulk")
    end

    push!(sim.interactions, dmi)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return dmi
end

""" 
    add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")

Add exchange interaction within the Geometry, or update corresponding A if other exch is added.
"""
function add_exch(sim::AbstractSim, geo::Geometry, A::Number; name="exch")
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.A, geo, A)
            return nothing
        end
    end
    n_total = sim.n_total
    field = zeros(Float64, 3 * n_total)
    energy = zeros(Float64, n_total)
    Spatial_A = zeros(Float64, n_total)
    update_scalar_geometry(Spatial_A, geo, A)
    if isa(sim, MicroSim)
        exch = Exchange(Spatial_A, field, energy, name)
    else
        exch = HeisenbergExchange(A, field, energy, name)
    end
    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.ite,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return exch
end

@doc raw"""
    add_exch_rkky(sim::AbstractSim, J::Float64; name="rkky")

Add an RKKY-type exchange to the system. The energy of RKKY-type exchange is defined as 

```math
E_\mathrm{rkky} =  - \int_\Gamma J_\mathrm{rkky} \mathbf{m}_{i} \cdot \mathbf{m}_{j} dA
```
where $\Gamma$ is the interface between two layers with magnetizations $\mathbf{m}_{i}$ and $\mathbf{m}_{j}$,
$J_\mathrm{rkky}$ is the coupling constant which is related to the spacer layer thickness. 

The effective field is given then as
```math
\mathbf{H}_i = \frac{1}{\mu_0 M_s}  \frac{J_\mathrm{rkky}}{\Delta} \mathbf{m}_{j} 
```
"""
function add_exch_rkky(sim::AbstractSim, J::Float64; name="rkky")
    n_total = sim.n_total
    field = zeros(Float64, 3 * n_total)
    energy = zeros(Float64, n_total)
    exch = ExchangeRKKY(J, field, energy, name)

    push!(sim.interactions, exch)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.ite,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return exch
end

"""
    add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)

Add Demag to the system. `Nx`, `Ny` and `Nz` can be used to describe the macro boundary conditions which means that
the given mesh is repeated `2Nx+1`, `2Ny+1 and `2Nz+1` times in `x`, `y` and `z` direction, respectively.
"""
function add_demag(sim::MicroSim; name="demag", Nx=0, Ny=0, Nz=0)
    demag = init_demag(sim, Nx, Ny, Nz)
    demag.name = name
    push!(sim.interactions, demag)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return demag
end

"""
    add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0,0,1), name="anis")

Add Anisotropy to the system, where the energy density is given by

```math
    E_\\mathrm{anis} =  K_{u} (1 - \\vec{m} \\cdot \\hat{u})^2
```
"""
function add_anis(sim::AbstractSim, Ku::NumberOrArrayOrFunction; axis=(0, 0, 1),
                  name="anis")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64
    Kus = zeros(T, n_total)
    init_scalar!(Kus, sim.mesh, Ku)

    Kus_kb = KernelAbstractions.zeros(backend[], T, n_total)
    field = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(backend[], T, n_total)

    lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)
    naxis = (axis[1] / lt, axis[2] / lt, axis[3] / lt)

    copyto!(Kus_kb, Kus)
    anis = Anisotropy(Kus_kb, naxis, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end

    @info "Uniaxial Anisotropy has been added."
    return anis
end

"""
    add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0,0,1), name="anis")

Add Anisotropy within the Geometry, or update corresponding Ku if other anis is added.
"""
#FIXME : fix this function
function add_anis(sim::AbstractSim, geo::Geometry, Ku::Number; axis=(0, 0, 1), name="anis")
    lt = sqrt(axis[1]^2 + axis[2]^2 + axis[3]^2)

    naxis = (axis[1] / lt, axis[2] / lt, axis[3] / lt)
    for interaction in sim.interactions
        if interaction.name == name
            update_scalar_geometry(interaction.Ku, geo, Ku)
            return nothing
        end
    end
    n_total = sim.n_total
    Kus = zeros(Float64, n_total)
    update_scalar_geometry(Kus, geo, Ku)
    field = zeros(Float64, 3 * n_total)
    energy = zeros(Float64, n_total)
    anis = Anisotropy(Kus, naxis, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end

"""
    update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name = "anis")

update anisotropy constant Ku according to its name.

Example:
```julia
    mesh = FDMesh(nx=200, ny=200, nz=12, dx=5e-9, dy=5e-9, dz=5e-9)
    sim = Sim(mesh)
    add_anis(sim, 3e4, axis = (0,0,1), name="K1")
    add_anis(sim, 1e5, axis = (1,0,0), name="K2")
    update_anis(sim, 5e4, name="K2")  #update anisotropy K2
```
"""
function update_anis(sim::MicroSim, Ku::NumberOrArrayOrFunction; name="anis")  #FIXME : fix this function
    n_total = sim.n_total
    Kus = zeros(Float64, n_total)
    init_scalar!(Kus, sim.mesh, Ku)
    for i in sim.interactions
        if i.name == name
            i.Ku[:] = Kus[:]
            return nothing
        end
    end
    return nothing
end

@doc raw"""
    add_cubic_anis(sim::AbstractSim, Kc::Float64; axis1=(1,0,0), axis2=(0,1,0), name="cubic")

add a cubic anisotropy with default axis (1,0,0) , (0,1,0), and (0,0,1). The third axis is defined as axis3 = axis1 x axis2.

```math
  E_\mathrm{cubic} = -\int_{V} K_c (m_x^4 + m_y^4 + m_z^4) \, dV
```

# Example
```julia
    add_cubic_anis(sim, 1e3, (1, 1, 0), (1, -1, 0))
```
"""
function add_cubic_anis(sim::AbstractSim, Kc::NumberOrArrayOrFunction; axis1=(1, 0, 0),
                        axis2=(0, 1, 0), name="cubic")
    n_total = sim.n_total
    T = single_precision.x ? Float32 : Float64
    Kcs = zeros(T, n_total)
    init_scalar!(Kcs, sim.mesh, Kc)

    norm1 = sqrt(axis1[1]^2 + axis1[2]^2 + axis1[3]^2)
    norm2 = sqrt(axis2[1]^2 + axis2[2]^2 + axis2[3]^2)
    naxis1, naxis2 = axis1 ./ norm1, axis2 ./ norm2
    if abs.(sum(naxis1 .* naxis2)) > 1e-10
        @error("The axis1 and axis2 are not perpendicular to each other!")
        return nothing
    end
    naxis3 = cross_product(axis1, axis2)

    Kcs_kb = KernelAbstractions.zeros(backend[], T, n_total)
    field = KernelAbstractions.zeros(backend[], T, 3 * n_total)
    energy = KernelAbstractions.zeros(backend[], T, n_total)

    copyto!(Kcs_kb, Kcs)

    anis = CubicAnisotropy(Kcs_kb, naxis1, naxis2, naxis3, field, energy, name)
    push!(sim.interactions, anis)

    if sim.save_data
        id = length(sim.interactions)
        push!(sim.saver.items,
              SaverItem(string("E_", name), "J",
                        o::AbstractSim -> sum(o.interactions[id].energy)))
    end
    return anis
end
