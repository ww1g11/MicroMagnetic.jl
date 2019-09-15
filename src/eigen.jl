using SparseArrays

function EigenProblem(mesh::Mesh, m0::TupleOrArrayOrFunction; J=1.0, mu_s=1, gamma=1, Kx=0, Kz=0)
    nxyz = mesh.nxyz
    m = zeros(Float64, 3*nxyz)
    init_vector!(m, mesh, m0)

    xyz = reshape(m, 3, nxyz)
    theta = zeros(Float64, nxyz)
    phi = zeros(Float64, nxyz)
    ct = zeros(Float64, nxyz) #cos(theta)
    st = zeros(Float64, nxyz) #sin(theta)
    cp = zeros(Float64, nxyz) #cos(phi)
    sp = zeros(Float64, nxyz) #sin(phi)
    for i = 1:nxyz
        r_xy = sqrt(xyz[1, i]^2 + xyz[2, i]^2)
        theta[i] = atan(r_xy, xyz[3, i])
        phi[i] = atan(xyz[2, i], xyz[1, i])
        ct[i] = cos(theta[i])
        st[i] = sin(theta[i])
        cp[i] = cos(phi[i])
        sp[i] = sin(phi[i])
        #println(abs(xyz[3, i]-ct[i]),"  ",abs(xyz[2, i]-st[i]*sp[i]), " ", abs(xyz[1, i]-st[i]*cp[i]))
    end
    A = spzeros(2*nxyz, 2*nxyz)

    bulid_matrix_J(A, J, mesh, phi, st, ct)
    bulid_matrix_K(A, Kx, Kz, mesh.nxyz, st, ct, sp, cp)

    return A
end

function bulid_matrix_J(A, J, mesh, phi, st, ct)
    n = mesh.nxyz


    for i = 1:n
        total = 0
        for k=1:6
            j = mesh.ngbs[k,i]
            if j<0
                continue;
            end

            total += J*(ct[i]*ct[j]+st[i]*st[j]*cos(phi[i]-phi[j]))

            #hu
            A[i+n, j] -= J*(cos(phi[i]-phi[j])*ct[i]*ct[j]+st[i]*st[j])
            A[i+n, j+n] -= J*ct[i]*sin(phi[i]-phi[j])

            #hv
            A[i, j] -= J*ct[j]*sin(phi[i]-phi[j])
            A[i, j+n] += J*cos(phi[i]-phi[j])
        end

        #H0
        A[i,i+n] -= total
        A[i+n,i] += total
    end
end

function bulid_matrix_K(A, Kx, Kz, nxyz, st, ct, sp, cp)
    n = nxyz

    for i = 1:n
        j = i + n
        #hu
        A[j,i] -= 2*Kx*(cp[i]*ct[i])^2 + 2*Kz*st[i]^2
        A[j,i] += 2*Kx*sp[i]*cp[i]*ct[i]

        #hv
        A[i,i] += -2*Kx*sp[i]*cp[i]*ct[i]
        A[i,j] += 2*Kx*sp[i]*sp[i]

        #H0
        total = 2*Kx*(cp[i]*st[i])^2+2*Kz*(ct[i])^2
        A[i,j] -= total
        A[j,i] += total
    end
end
