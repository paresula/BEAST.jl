struct HH3DSphericalWave{T,P} <: Functional
    location::P
    wavenumber::T
    amplitude::T
end

function HH3DSphericalWave(direction, wavenumber; amplitude=1)
    l, a = promote(wavenumber, amplitude)
    HH3DSphericalWave(direction, l, a)
end

function (f::HH3DSphericalWave)(r)
    l = f.location 
    k = f.wavenumber
    a = f.amplitude
    a * exp(-im*k*norm(r-l))/(4*π*norm(r-l))
end

struct HH3DPlaneWave{T,P} <: Functional
    direction::P
    wavenumber::T
    amplitude::T
end

function HH3DPlaneWave(direction, wavenumber; amplitude=1)
    w, a = promote(wavenumber, amplitude)
    HH3DPlaneWave(direction, w, a)
end

function (f::HH3DPlaneWave)(r)
    d = f.direction
    k = f.wavenumber
    a = f.amplitude
    a * exp(-im*k*dot(d,r))
end

struct NormalDerivative{F} <: Functional
    field::F
end

const ∂n = Val{:normalderivative}
(::Type{Val{:normalderivative}})(f) = NormalDerivative(f)

function (f::NormalDerivative{T})(manipoint) where T<:HH3DPlaneWave
    d = f.field.direction
    k = f.field.wavenumber
    a = f.field.amplitude
    n = normal(manipoint)
    r = cartesian(manipoint)
    -im*k*a * dot(d,n) * exp(-im*k*dot(d,r))
end

function (f::NormalDerivative{T})(manipoint) where T<:HH3DSphericalWave
    l = f.field.location
    k = f.field.wavenumber
    a = f.field.amplitude
    n = normal(manipoint)
    r0 = cartesian(manipoint)
    r = norm(l-r0)
    er = (l-r0)/r
    dot(exp(-im*k*r)*(-1-i*k*r)/(4*π*r^2)*er,n)
end

integrand(::NormalDerivative, test_vals, field_vals) = dot(test_vals[1], field_vals)


#= abstract type Charge end

mutable struct HH3DPointCharge{T} <: Charge
    location::SVector{3,T}
    charge::T
end

function HH3DPointCharge(l, c)
    T = promote_type(eltype(l), typeof(c))
    HH3DPointCharge{T}(l, c)
end =#