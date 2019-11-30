import Base.+
import Base.-
import Base.*

struct Epsilon
    x0::Float64  #the constant
    x1::Float64  #the coefficient of ui
    x2::Float64  #the coefficient of vi
    x3::Float64  #the coefficient of uj
    x4::Float64  #the coefficient of vj
end

struct VecEps
    x::Epsilon
    y::Epsilon
    z::Epsilon
end

function +(v1::Epsilon, v2::Epsilon)
    x0 = v1.x0 + v2.x0
    x1 = v1.x1 + v2.x1
    x2 = v1.x2 + v2.x2
    x3 = v1.x3 + v2.x3
    x4 = v1.x4 + v2.x4
    return Epsilon(x0, x1, x2, x3, x4)
end

function -(v1::Epsilon, v2::Epsilon)
    x0 = v1.x0 - v2.x0
    x1 = v1.x1 - v2.x1
    x2 = v1.x2 - v2.x2
    x3 = v1.x3 - v2.x3
    x4 = v1.x4 - v2.x4
    return Epsilon(x0, x1, x2, x3, x4)
end

function *(v1::Epsilon, v2::Epsilon)
    x0 = v1.x0 * v2.x0
    x1 = v1.x0 * v2.x1 + v2.x0 * v1.x1
    x2 = v1.x0 * v2.x2 + v2.x0 * v1.x2
    x3 = v1.x0 * v2.x3 + v2.x0 * v1.x3
    x4 = v1.x0 * v2.x4 + v2.x0 * v1.x4
    return Epsilon(x0, x1, x2, x3, x4)
end

function *(v1::Number, v2::Epsilon)
    x0 = v1 * v2.x0
    x1 = v1 * v2.x1
    x2 = v1 * v2.x2
    x3 = v1 * v2.x3
    x4 = v1 * v2.x4
    return Epsilon(x0, x1, x2, x3, x4)
end

function eps_cross(v1::VecEps, v2::VecEps)
    x =  v1.y*v2.z - v1.z*v2.y
    y =  v1.z*v2.x - v1.x*v2.z
    z =  v1.x*v2.y - v1.y*v2.x
    return VecEps(x, y, z)
end


function test_epsilon()
    mxi = Epsilon(0,1.0,0,0,0)  #local spin i
    myi = Epsilon(0,0,1.0,0,0)
    mzi = Epsilon(1.0,0,0,0,0)
    mxj = Epsilon(0,0,0,1.0,0)  #local spin j
    myj = Epsilon(0,0,0,0,1.0)
    mzj = Epsilon(1.0,0,0,0,0)
    mi = VecEps(mxi, myi, mzi)
    mj = VecEps(mxj, myj, mzj)
    println(cross(mi, mj))
end
