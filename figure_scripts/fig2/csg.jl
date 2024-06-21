
using JuMag

mesh = FDMesh(dx=2e-9, dy=2e-9, dz=2e-9, nx=100, ny=100, nz=50)

t1 = Torus(R = 60e-9, r=20e-9)
save_vtk(mesh, t1, "shape1")

b1 = Box(sides = (160e-9, 50e-9, 40e-9), theta=pi/4)
t2 = t1 + b1
save_vtk(mesh, t2, "shape2")

t3 = t1 - b1 
save_vtk(mesh, t3, "shape3")

t4 = t1 * b1 
save_vtk(mesh, t4, "shape4")





