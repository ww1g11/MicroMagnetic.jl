using JuMag
using Test
#using LinearAlgebra
using Arpack

function energy(n, J=1, K=0.01, L=100)
    return 2*K + J*(n*pi)^2/L^2
end

mesh =  CubicMesh(nx=100, ny=1, nz=1, a=1.0, pbc="open")
A = EigenProblem(mesh, (1,0,0), J=1.0, Kx=0.01)
E, x = eigs(A, nev=20, which=:LM, sigma=1e-8)
E = imag(E)
Es = filter((x) -> x>0, E)
println(Es)
@test isapprox(Es[1],energy(0))
@test isapprox(Es[2],energy(1), atol=1e-5)
@test isapprox(Es[3],energy(2), atol=1e-5)
