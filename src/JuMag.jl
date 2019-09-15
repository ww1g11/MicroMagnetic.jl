__precompile__()

module JuMag
using Printf

export init_m0,
       add_zeeman,
       add_dmi,
       add_exch, add_anis, add_demag, add_exch_rkky,
       update_zeeman,update_anis,add_exch_vector,
       run_until, relax,
       save_vtk, save_m,
       FDMesh, set_Ms, Sim,
       CubicMesh, set_mu_s,
       set_ux, write_data,
       compute_system_energy,
       compute_skyrmion_number,
       compute_guiding_centre, set_aj,
       Neb,interpolate_m,save_ovf,read_ovf,
       EigenProblem

export mu_0, mu_B, k_B, c_e, eV, meV, m_e, g_e, h_bar, gamma, mu_s_1, h_bar_gamma, mT

const _cuda_using_double = Ref(false)
const _cuda_available = Ref(true)
const _using_gpu = Ref(false)

function cuda_using_double(flag = true)
   _cuda_using_double[] = flag
   return nothing
end

include("const.jl")
include("head.jl")
include("util.jl")
include("mesh.jl")
include("driver.jl")
include("sd.jl")
include("llg.jl")
include("rk.jl")
include("helper.jl")
include("ode.jl")
include("fileio.jl")
include("sim.jl")
include("demag.jl")
include("vtk.jl")
include("neb.jl")
include("eigen.jl")
include("ovf2.jl")

try
	using CUDAnative, CuArrays
    #CuArrays.allowscalar(false) TODO: where should it be?
	#using CuArrays.CUFFT
	#@info "Running CUFFT $(CUFFT.version())"
catch
    _cuda_available[] = false
    _using_gpu[] = false
    @warn "CUDA is not available!"
end

if _cuda_available.x
    include("cuda/head.jl")
    include("cuda/mesh.jl")
    include("cuda/driver.jl")
    include("cuda/sim.jl")
    include("cuda/ode.jl")
    include("cuda/llg.jl")
    include("cuda/util.jl")
    include("cuda/kernels.jl")
    include("cuda/field.jl")
    include("cuda/demag_kernel.jl")
    include("cuda/demag.jl")
    include("cuda/sd.jl")
    include("cuda/mc.jl")
    include("cuda/mc_kernel.jl")
    include("cuda/vtk.jl")
    include("cuda/ovf2.jl")
    export FDMeshGPU, CubicMeshGPU, TriangularMesh, MonteCarlo, run_sim, add_demag_gpu
end


end
