abstract type NEBDriver end
mutable struct NEBSaver
  name::String
  nsteps::Int64
  header_saved::Bool
  headers::Array  #string or tuple<string> array
  units::Array #string or tuple<string> array
  results::Array  #function array
end
mutable struct NEB
  N::Int64
  sim::AbstractSim
  init_m::Any
  intervals::Any
  driver::NEBDriver
  images::Array{Float64, 2}
  pre_images::Array{Float64, 2}
  name::String
  field::Array{Float64, 2}
  energy::Array{Float64, 2}
  saver::NEBSaver
end
function Neb(sim::AbstractSim,init_m::Any, intervals::Any;name="NEB", driver = "NEB_SD")
  nxyz = sim.mesh.nx*sim.mesh.ny*sim.mesh.nz
  N = sum(intervals)+length(intervals)+1
  images = zeros(Float64,3*sim.nxyz,N)
  pre_images = zeros(Float64,3*sim.nxyz,N)
  energy = zeros(Float64,nxyz,N)
  Ms = zeros(Float64,nxyz)
  driver = create_neb_driver(driver, nxyz,N)
  field = zeros(Float64,3*nxyz,N)
  energy = zeros(nxyz,N)

  headers = ["steps"]
  units = ["<>"]
  results = Any[o::NEB -> o.driver.steps]
  saver = NEBSaver(string(name, ".txt"), 0, false, headers, units, results)
  #interactions = []
  neb = NEB(N,sim,init_m,intervals,driver,images,pre_images,name,field,energy,saver)
  init_images(neb)
  return neb
end
function create_neb_driver(driver::String, nxyz::Int64,N::Int64) #TODO: FIX ME
  if driver=="NEB_SD"
      gk = zeros(Float64,3*nxyz,N)
  return NEB_SD(gk, 0.0, 1e-4, 1e-14, 0)
  else
    println("NEB_SD supported!")
  end
end
mutable struct NEB_SD <: NEBDriver
  gk::Array{Float64, 2}
  tau::Float64
  max_tau::Float64
  min_tau::Float64
  steps::Int64
end
@inline function cross_product(a::Array{T,1}, b::Array{T,1}) where {T<:AbstractFloat}
  x1=a[1];x2=a[2];x3=a[3];y1=b[1];y2=b[2];y3=b[3]
  return [-x3*y2 + x2*y3, x3*y1 - x1*y3, -x2*y1 + x1*y2]
end
function rotating_oprate(m1::Array{Float64,1},m2::Array{Float64,1},theta::Float64) ##return m'=retate(m1,theta)
  m=cross_product(m1,m2)
  if m==[0,0,0]
    if m1==[0,0,1.0]
      m= [1.0,0,0]
    else
     m=cross_product(m1,[0,0,1.0])
    end
  end
  normalise(m,1)
  op=[cos(theta)+(1-cos(theta))*m[1]^2 (1-cos(theta))*m[1]*m[2]-sin(theta)*m[3] (1-cos(theta))*m[1]*m[3]+sin(theta)*m[2];(1-cos(theta))*m[2]*m[3]+sin(theta)*m[3] cos(theta)+(1-cos(theta))*m[2]*m[2] (1-cos(theta))*m[2]*m[3]-sin(theta)*m[1];(1-cos(theta))*m[1]*m[3]-sin(theta)*m[2] (1-cos(theta))*m[2]*m[3]+sin(theta)*m[1] cos(theta)+(1-cos(theta)*m[3]*m[3])]
  return op*m1
end
function interpolate_m(m1::Array{Float64,1}, m2::Array{Float64,1}, n::Int)
  nxyz=Int(length(m1)/3)
  normalise(m1,nxyz)
  normalise(m2,nxyz)
  m = zeros(3*nxyz,n+1)
  b1=reshape(m1,3,nxyz)
  b2=reshape(m2,3,nxyz)
  for i=1:n+1
    for j=1:nxyz
      k=3*j-2
      theta=acos(b1[:,j]'*b2[:,j])
      dtheta=theta/(n+1)
      angle=(i-1)dtheta
      m[k,i],m[k+1,i],m[k+2,i]=rotating_oprate(b1[:,j],b2[:,j],angle)
    end
  end
  return m
end
function init_images(neb::NEB)
  sim=neb.sim
  nxyz=sim.nxyz
  N=neb.N
  m=zeros(3*nxyz,N)
  pics= neb.init_m
  intervals=neb.intervals
  if length(pics)==length(intervals)+1
    n=1
    for i=1:(length(pics)-1)
      m1=zeros(3*nxyz)
      m2=zeros(3*nxyz)
      init_vector!(m1, sim.mesh, pics[i])
      init_vector!(m2, sim.mesh, pics[i+1])
      M=interpolate_m(m1,m2,Int(intervals[i]))
      for j=1:intervals[i]+1
          m[:,n]=M[:,j]
          n+=1
      end
    end
    m3=zeros(3*nxyz)
    init_vector!(m3, sim.mesh, pics[length(pics)])
    m[:,N]=m3[:]
  else
    println("Input error!")
  end
  neb.pre_images[:]=m[:]
  neb.images[:]=neb.pre_images[:]
  for n = 1:N
    name = @sprintf("E_total_%g",n)
    push!(neb.saver.headers,name)
    push!(neb.saver.units, "J")
    fun =  o::NEB -> sum(o.energy[:,n])
    push!(neb.saver.results, fun)
  end
end
function compute_tau(driver::NEB_SD, pre_images::Array{Float64, 2}, images::Array{Float64, 2}, h::Array{Float64, 2}, nxyz::Int64,N)

  if driver.steps == 0
    for n=1:N, i=0:nxyz-1  #compute gk[] for step 0
        j = 3*i+1
        fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
        gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
        driver.gk[j,n] = gx
        driver.gk[j+1,n] = gy
        driver.gk[j+2,n] = gz
      end
     driver.tau  = driver.min_tau
     return nothing
  end
  sum1, sum2, sum3 = zeros(N),zeros(N),zeros(N)
  for n=1:N
    for i=0:nxyz-1
        j = 3*i+1
        sx = images[j,n] - pre_images[j,n]
        sy = images[j+1,n] - pre_images[j+1,n]
        sz = images[j+2,n] - pre_images[j+2,n]
        fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
        gx,gy,gz = cross_product(images[j,n],images[j+1,n],images[j+2,n], fx,fy,fz)
        fx = gx - driver.gk[j,n]
        fy = gy - driver.gk[j+1,n]
        fz = gz - driver.gk[j+2,n]
        driver.gk[j,n] = gx
        driver.gk[j+1,n] = gy
        driver.gk[j+2,n] = gz
        sum1[n] += sx*sx+sy*sy+sz*sz
        sum2[n] += sx*fx+sy*fy+sz*fz
        sum3[n] += fx*fx+fy*fy+fz*fz
     end
  end
     sum_1,sum_2,sum_3 = sum(sum1)/N,sum(sum2)/N,sum(sum3)/N
     tau1 = sum_2!=0.0 ? sum_1/sum_2 : driver.min_tau
     tau2 = sum_3!=0.0 ? sum_2/sum_3 : driver.min_tau
	   driver.tau = driver.steps%2 == 0 ? abs(tau2) : abs(tau1)
     if driver.tau > driver.max_tau
        driver.tau = driver.max_tau
     end
     return nothing
end
function run_step(neb::NEB)
  driver = neb.driver
  N = neb.N
  images = neb.images
  pre_images = neb.pre_images
  sim = neb.sim
  nxyz = sim.nxyz
  effective_field(neb, 0.0)
  compute_tau(driver, pre_images, images, neb.field, nxyz, N)
  neb.pre_images[:] =  neb.images[:]
  h = neb.field
  gk = driver.gk
  tau = driver.tau
  for n =1:N
    if (n == 1)||(n==N)
      continue
    end
    for i=0:sim.nxyz-1
      j = 3*i+1
      fx,fy,fz = cross_product(images[j,n],images[j+1,n],images[j+2,n], h[j,n],h[j+1,n],h[j+2,n])
      factor = 0.25*(fx*fx+fy*fy+fz*fz)*tau^2
      mx = (1-factor)*images[j,n] - tau*gk[j,n]
      my = (1-factor)*images[j+1,n] - tau*gk[j+1,n]
      mz = (1-factor)*images[j+2,n] - tau*gk[j+2,n]
      images[j,n] = mx/(1+factor)
      images[j+1,n] = my/(1+factor)
      images[j+2,n] = mz/(1+factor)
    end
  end
  for n=1:N
  max_length_error = error_length_m(neb.images[:,n], sim.nxyz)
  if max_length_error > 1e-15
    neb.images[:,n] = normalized(neb.images[:,n], sim.nxyz)                                 #fix me
  end
  end
  driver.steps += 1
  return  nothing
end
function normalized(a::Array{T, 1}, N::Int64) where {T<:AbstractFloat}
  for i = 0:N-1
      j = 3*i+1
      length = sqrt(a[j]*a[j]+a[j+1]*a[j+1]+a[j+2]*a[j+2])
      if length > 0
        length = 1.0/length
        a[j] *= length
        a[j+1] *= length
        a[j+2] *= length
      end
    end
   return a
end
function effective_field(neb::NEB,  t::Float64)
  sim = neb.sim
  nxyz = sim.nxyz
  images = neb.images
  N = neb.N
  fill!(neb.field, 0.0)
  fill!(neb.energy, 0.0)
  for n = 1:neb.N
  for interaction in sim.interactions
    effective_field(interaction, sim, images[:,n], t)
    neb.field[:,n] .+= interaction.field[:]
    neb.energy[:,n] .+= interaction.energy[:]
  end
  end
  t = zeros(3*nxyz,N)
  for n = 1:neb.N
    if (n==1)||(n==neb.N)
      continue
    end
    E1 = sum(neb.energy[:,n-1])
    E2 = sum(neb.energy[:,n])
    E3 = sum(neb.energy[:,n+1])
    if (E1>E2)&&(E2>E3)
      t[:,n] = normalized((images[:,n-1]-images[:,n]),nxyz)
    elseif (E3>E2)&&(E2>E1)
      t[:,n] = normalized(images[:,n+1]-images[:,n],nxyz)
    else
      t[:,n] = normalized(images[:,n-1]-images[:,n],nxyz)
    end
    neb.field[:,n] = neb.field[:,n] - neb.field[:,n]*t[:,n]'*t[:,n]
  end
  return 0
end
function compute_system_energy(neb::NEB,  t::Float64)
  sim = neb.sim
  images = neb.images
  #sim.total_energy = 0
  fill!(neb.energy, 0.0)
  for n = 1:neb.N
  for interaction in sim.interactions
    effective_field(interaction, sim, images[:,n], t)
	neb.energy[:,n] .+= interaction.energy[:]
  end
  end
  return 0
end
function relax(neb::NEB; maxsteps=10000, stopping_dmdt=0.01, stopping_torque=0.1, save_m_every = 10, save_vtk_every=-1, vtk_folder="vtks")

  if isa(neb.driver, NEB_SD)
      is_relax_NEB = true
  end

  if !isdir(vtk_folder)
      mkdir(vtk_folder)
  end

  if is_relax_NEB
      relax_NEB(neb, maxsteps, Float64(stopping_torque), save_m_every, save_vtk_every, vtk_folder)
  end
  return nothing
end
function relax_NEB(neb::NEB, maxsteps::Int64, stopping_torque::Float64, save_m_every::Int64, save_vtk_every::Int64, vtk_folder::String)
  N = neb.N
  sim = neb.sim
  gk_abs = zeros(Float64,3*sim.nxyz)
  driver = neb.driver
  maxtorque=zeros(N)
  if save_m_every>0  
    compute_system_energy(neb,  0.0)
    write_data(neb)
  end
  if save_vtk_every > 0 
    save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, 0)))
  end
  for i=1:maxsteps
      run_step(neb)
      #maxtorque = zeros(N)
      for n = 1:N
        if (n==1)||(n==N)
          continue
        end
        abs!(gk_abs, driver.gk[:,n])  #max_torque = maximum(abs.(driver.gk)) eats gpu memory???
        maxtorque[n] = maximum(gk_abs)
      end
      max_torque = maximum(maxtorque)
      @info @sprintf("step=%5d  tau=%10.6e  max_torque=%10.6e", i, driver.tau, max_torque)
      if save_m_every>0 && i%save_m_every == 0
        compute_system_energy(neb,  0.0)
        write_data(neb)
      end
      if save_vtk_every > 0 && i%save_vtk_every == 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
      end
      if max_torque < stopping_torque
      @info @sprintf("max_torque (mxmxH) is less than stopping_torque=%g, Done!", stopping_torque)
      if save_m_every>0
        compute_system_energy(neb, 0.0)
        write_data(neb)
      end
      if save_vtk_every > 0
        save_vtk(neb, joinpath(vtk_folder, @sprintf("%s_%d", neb.name, i)))
      end
      break
      end
  end
return nothing
end
function write_data(neb::NEB)
  sim = neb.sim
  saver = neb.saver
  if !saver.header_saved
    name = @sprintf("%s", saver.name)
    io = open(name, "w");
    write(io, "#")
    for header in saver.headers
      write(io, formatstring(header));
    end
    write(io, "\n#");

    for s in saver.units
      write(io, formatstring(s));
    end
    write(io, "\n")
    saver.header_saved = true
  else
    io = open(saver.name, "a");
  end
  for fun in saver.results
    write(io, formatstring(fun(neb)));
  end
  write(io, "\n")
  close(io);
end
function save_vtk(neb::NEB, fname::String; fields::Array{String, 1} = String[])
  sim = neb.sim
  mesh = sim.mesh
  nx, ny, nz = mesh.nx, mesh.ny, mesh.nz
  xyz = zeros(Float32, 3, nx+1, ny+1, nz+1)
  dx, dy, dz = 1.0,1.0,1.0
  if isa(mesh, CubicMesh)
    dx, dy, dz=mesh.a, mesh.a, mesh.a
  else
    dx, dy, dz=mesh.dx, mesh.dy, mesh.dz
  end
  for k = 1:nz+1, j = 1:ny+1, i = 1:nx+1
    xyz[1, i, j, k] = (i-0.5)*dx
    xyz[2, i, j, k] = (j-0.5)*dy
    xyz[3, i, j, k] = (k-0.5)*dz
  end
  for n =1:neb.N
    name = @sprintf("%s_%d", fname, n)
    vtk = vtk_grid(name, xyz)
  b = reshape(neb.images[:,n], (3, nx, ny, nz))
  vtk_cell_data(vtk, b , "m")
  if length(fields) > 0
    fields = Set(sim.fields)
    for i in sim.interactions
      if i.name in fields
        effective_field(i, sim, neb.images[:,n], 0.0)
        b = reshape(i.field, (3, nx, ny, nz))
        vtk_cell_data(vtk, b, i.name)
      end
    end
  end
  vtk_save(vtk)
  end
end
