
struct HH3DDoubleLayerNear{K}
  gamma::K
end

struct HH3DSingleLayerNear{K}
  gamma::K
end

function HH3DDoubleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DDoubleLayerNear(-imag(wavenumber))
  else
    HH3DDoubleLayerNear(wavenumber*im)
  end
end

function HH3DSingleLayerNear(;wavenumber=error("wavenumber is a required argument"))
  if iszero(real(wavenumber))
    HH3DSingleLayerNear(-imag(wavenumber))
  else
    HH3DSingleLayerNear(wavenumber*im)
  end
end

const HHNearField3D = Union{HH3DDoubleLayerNear, HH3DSingleLayerNear}
quaddata(op::HHNearField3D,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::HHNearField3D,refspace,p,y,q,el,qdata) = qdata[1,q]

#= quaddata(op::HH3DSingleLayerNear,rs,els) = quadpoints(rs,els,(3,))
quadrule(op::HH3DSingleLayerNear,refspace,p,y,q,el,qdata) = qdata[1,q] 
 =#
function kernelvals(op::HHNearField3D,y,p)

  γ = op.gamma
  x = cartesian(p)
  r = y - x
  R = norm(r)
  γR = γ*R

  inv_R = 1/R

  expn = exp(-γR)
  green = expn * inv_R * inv_4pi
  gradgreen = -(γ + inv_R) * green * inv_R * r

  nx = normal(p)

  (;γ, r, R, green, gradgreen, nx)
end

function integrand(op::HH3DDoubleLayerNear,krn,y,f,p)

  ∇G = krn.gradgreen
  nx = krn.nx

  fx = f.value

  ∂G∂n = nx ⋅ ∇G

  return ∂G∂n * fx 
end

#= function kernelvals(op::HH3DSingleLayerNear, y, p)

  γ = op.gamma
  x = cartesian(p)
  r = y - x
  r = y - x
  R = norm(r)
  γR = γ*R
  inv_R = 1/R
  expn = exp(-γR)
  green = expn * inv_R *inv_4pi
  (;γ, r, R, green)

end =#

function integrand(op::HH3DSingleLayerNear, krn, y, f, p)
  G = krn.green
  fx = f.value
  return G*fx
end