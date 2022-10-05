using BEAST
using CompScienceMeshes
using StaticArrays
using LinearAlgebra
using Test
using Plotly

## Looking at convergence 
hs = [0.5,0.45,0.4,0.35,0.3,0.25,0.2,0.15,0.1,0.08, 0.07]
#irs = [0.99,0.95,0.91,0.87, 0.83, 0.79, 0.75]
irs=[0.8]
err_IDPSL_pot=zeros(Float64, length(hs), length(irs))
err_IDPDL_pot=zeros(Float64, length(hs), length(irs))
err_INPSL_pot=zeros(Float64, length(hs), length(irs))
err_INPDL_pot=zeros(Float64, length(hs), length(irs))
err_IDPSL_field=zeros(Float64, length(hs), length(irs))
err_IDPDL_field=zeros(Float64, length(hs), length(irs))
err_INPSL_field=zeros(Float64, length(hs), length(irs))
err_INPDL_field=zeros(Float64, length(hs), length(irs))

err_EDPSL_pot=zeros(Float64, length(hs), length(irs))
err_EDPDL_pot=zeros(Float64, length(hs), length(irs))
err_ENPSL_pot=zeros(Float64, length(hs), length(irs))
err_ENPDL_pot=zeros(Float64, length(hs), length(irs))

err_EDPSL_field=zeros(Float64, length(hs), length(irs))
err_EDPDL_field=zeros(Float64, length(hs), length(irs))
err_ENPSL_field=zeros(Float64, length(hs), length(irs))
err_ENPDL_field=zeros(Float64, length(hs), length(irs))
for (i,h) in enumerate(hs)

    r=50.0
    sphere = meshsphere(r,h*r)
    X0 = lagrangecxd0(sphere)
    X1 = lagrangec0d1(sphere)

    𝓢 = Helmholtz3D.singlelayer(gamma=0.0)
    𝓓 = Helmholtz3D.doublelayer(gamma=0.0)
    𝓓t = Helmholtz3D.doublelayer_transposed(gamma=0.0)
    𝓝 = -Helmholtz3D.hypersingular(gamma=0.0)

    q=100.0
    ϵ=1.0

    pos1 = SVector(r*1.5,0.0,0.0)
    pos2 = SVector(-r*1.5,0.0,0.0)
    Φ_inc(x) = q/(4*π*ϵ)*(1/(norm(x-pos1))-1/(norm(x-pos2)))
    ∂nΦ_inc(x) = -q/(r*4*π*ϵ)*((norm(x)^2-dot(pos1,x))/(norm(x-pos1)^3)-(norm(x)^2-dot(pos2,x))/(norm(x-pos2)^3))
    fieldtheo(x) = q/(4*π*ϵ)*((x-pos1)/(norm(x-pos1)^3)-(x-pos2)/(norm(x-pos2)^3))

    gD0 = assemble(ScalarTrace(Φ_inc),X0)
    gD1 = assemble(ScalarTrace(Φ_inc),X1)
    gN = assemble(ScalarTrace(∂nΦ_inc), X1)

    G = assemble(Identity(), X1, X1)
    𝗼=ones(numfunctions(X1))

    M_IDPSL = assemble(𝓢, X0, X0)
    M_IDPDL = (-1/2*assemble(Identity(),X1,X1) + assemble(𝓓, X1,X1))

    M_INPDL = -assemble(𝓝, X1, X1)+G*𝗼*𝗼'*G
    M_INPSL = (1/2*assemble(Identity(),X1,X1) + assemble(𝓓t, X1, X1))

    ρ_IDPSL = M_IDPSL \ (-gD0)
    ρ_IDPDL = M_IDPDL \ (gD1)

    ρ_INPDL = M_INPDL \ (-gN)
    ρ_INPSL = M_INPSL \ (-gN)

    #plot(patch(sphere, real.(facecurrents(ρ1,X0)[1])))
    for (j,ir) in enumerate(irs)
    pts = meshsphere(r*ir, r*ir*0.6).vertices

    pot_IDPSL = potential(HH3DSingleLayerNear(0.0), pts,ρ_IDPSL, X0, type=ComplexF64)
    pot_IDPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_IDPDL, X1, type = ComplexF64)

    pot_INPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_INPSL, X1, type=ComplexF64)
    pot_INPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_INPDL, X1, type = ComplexF64)

    err_IDPSL_pot[i,j] = norm(pot_IDPSL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_IDPDL_pot[i,j] = norm(pot_IDPDL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_INPSL_pot[i,j] = norm(pot_INPSL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_INPDL_pot[i,j] = norm(pot_INPDL+Φ_inc.(pts))./norm(Φ_inc.(pts))

    field_IDPSL = potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_IDPSL, X0)
    field_IDPDL = potential(HH3DHyperSingularNear(0.0), pts, ρ_IDPDL, X1)
    field_INPSL = potential(HH3DDoubleLayerTransposedNear(0.0), pts, ρ_INPSL, X1)
    field_INPDL = potential(HH3DHyperSingularNear(0.0), pts, ρ_INPDL, X1)
  
    err_IDPSL_field[i,j] = norm(field_IDPSL-fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_IDPDL_field[i,j] = norm(field_IDPDL+fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_INPSL_field[i,j] = norm(field_INPSL-fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_INPDL_field[i,j] = norm(field_INPDL+fieldtheo.(pts))/norm(fieldtheo.(pts))
    

    end
    pos1 = SVector(r*0.5,0.0,0.0)
    pos2 = SVector(-r*0.5,0.0,0.0)
    Φ_inc(x) = q/(4*π*ϵ)*(1/(norm(x-pos1))-1/(norm(x-pos2)))
    ∂nΦ_inc(x) = -q/(r*4*π*ϵ)*((norm(x)^2-dot(pos1,x))/(norm(x-pos1)^3)-(norm(x)^2-dot(pos2,x))/(norm(x-pos2)^3))
    fieldtheo(x) = q/(4*π*ϵ)*((x-pos1)/(norm(x-pos1)^3)-(x-pos2)/(norm(x-pos2)^3))

    gD0 = assemble(ScalarTrace(Φ_inc),X0)
    gD1 = assemble(ScalarTrace(Φ_inc),X1)
    gN = assemble(ScalarTrace(∂nΦ_inc), X1)

    G = assemble(Identity(), X1, X1)
    𝗼=ones(numfunctions(X1))

    M_EDPSL = assemble(𝓢, X0, X0)
    M_EDPDL = (1/2*assemble(Identity(),X1,X1) + assemble(𝓓, X1,X1))

    M_ENPDL = -assemble(𝓝, X1, X1)+G*𝗼*𝗼'*G
    M_ENPSL = (-1/2*assemble(Identity(),X1,X1) + assemble(𝓓t, X1, X1))

    ρ_EDPSL = M_EDPSL \ (-gD0)
    ρ_EDPDL = M_EDPDL \ (gD1)

    ρ_ENPDL = M_ENPDL \ (-gN)
    ρ_ENPSL = M_ENPSL \ (-gN)

    #plot(patch(sphere, real.(facecurrents(ρ1,X0)[1])))
    for (j,ir) in enumerate(irs)
    testsphere = meshsphere(r/ir, r/ir*0.6)
    pts = testsphere.vertices[norm.(testsphere.vertices).>r]

    pot_EDPSL = potential(HH3DSingleLayerNear(0.0), pts,ρ_EDPSL, X0, type=ComplexF64)
    pot_EDPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_EDPDL, X1, type = ComplexF64)

    pot_ENPSL = potential(HH3DSingleLayerNear(0.0), pts, ρ_ENPSL, X1, type=ComplexF64)
    pot_ENPDL = potential(HH3DDoubleLayerNear(0.0), pts, ρ_ENPDL, X1, type = ComplexF64)

    err_EDPSL_pot[i,j] = norm(pot_EDPSL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_EDPDL_pot[i,j] = norm(pot_EDPDL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_ENPSL_pot[i,j] = norm(pot_ENPSL+Φ_inc.(pts))./norm(Φ_inc.(pts))
    err_ENPDL_pot[i,j] = norm(pot_ENPDL+Φ_inc.(pts))./norm(Φ_inc.(pts))

    field_EDPSL = -potential(HH3DDoubleLayerTransposedNear(im*k), pts, ρ_EDPSL, X0)
    field_EDPDL = potential(HH3DHyperSingularNear(im*k), pts, ρ_EDPDL, X1)
    field_ENPSL = -potential(HH3DDoubleLayerTransposedNear(im*k), pts, ρ_ENPSL, X1)
    field_ENPDL = potential(HH3DHyperSingularNear(im*k), pts, ρ_ENPDL, X1)

    err_EDPSL_field[i,j] = norm(field_EDPSL+fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_EDPDL_field[i,j] = norm(field_EDPDL+fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_ENPSL_field[i,j] = norm(field_ENPSL+fieldtheo.(pts))/norm(fieldtheo.(pts))
    err_ENPDL_field[i,j] = norm(field_ENPDL+fieldtheo.(pts))/norm(fieldtheo.(pts))
end
end

##

plot([  scatter(x = hs, y = err_IDPSL_pot[:,end], name = "IDPSL"),
        scatter(x = hs, y = err_IDPDL_pot[:,end], name = "IDPDL"),
        scatter(x = hs, y = err_INPSL_pot[:,end], name = "INPSL"),
        scatter(x = hs, y = err_INPDL_pot[:,end], name = "INPDL"),
        scatter(x=hs,y=hs.^2/1000, name = "h^2"),
        scatter(x=hs,y=hs.^3/10, name = "h^3"),
        scatter(x=hs,y=hs, name = "h")],
        Layout(
            yaxis_type="log", 
            xaxis_type="log",
            title = "Errors of the potential for the interior problem formualtions")
    )

plot([  scatter(x = hs, y = err_EDPSL_pot[:,end], name = "EDPSL"),
    scatter(x = hs, y = err_EDPDL_pot[:,end], name = "EDPDL"),
    scatter(x = hs, y = err_ENPSL_pot[:,end], name = "ENPSL"),
    scatter(x = hs, y = err_ENPDL_pot[:,end], name = "ENPDL"),
    scatter(x=hs,y=hs.^2/hs[end]^2*err_ENPDL_pot[end], name = "h^2"),
    #scatter(x=hs,y=hs.^3/10, name = "h^3"),
    scatter(x=hs,y=hs/hs[end]*err_EDPDL_pot[end], name = "h")],
    Layout(
        yaxis_type="log", 
        xaxis_type="log",
        title = "Errors of the potential for the exterior problem formualtions")
)

plot([  scatter(x = hs, y = err_IDPSL_field[:,end], name = "IDPSL"),
        scatter(x = hs, y = err_IDPDL_field[:,end], name = "IDPDL"),
        scatter(x = hs, y = err_INPSL_field[:,end], name = "INPSL"),
        scatter(x = hs, y = err_INPDL_field[:,end], name = "INPDL"),
        scatter(x=hs,y=hs.^2/hs[end]^2*err_INPDL_field[end], name = "h^2"),
        scatter(x=hs,y=hs.^3/hs[end]^3*err_IDPSL_field[end], name = "h^3"),
        scatter(x=hs,y=hs/hs[end]*err_IDPDL_field[end], name = "h")],
        Layout(
            yaxis_type="log", 
            xaxis_type="log",
            title = "Errors of the field for the interior problem formualtions")
    )

plot([  scatter(x = hs, y = err_EDPSL_field[:,end], name = "EDPSL"),
    scatter(x = hs, y = err_EDPDL_field[:,end], name = "EDPDL"),
    scatter(x = hs, y = err_ENPSL_field[:,end], name = "ENPSL"),
    scatter(x = hs, y = err_ENPDL_field[:,end], name = "ENPDL")],
    #scatter(x=hs,y=hs.^2/hs[end]^2*err_ENPDL_field[end], name = "h^2"),
    #scatter(x=hs,y=hs.^3/hs[end]^3*err_EDPSL_field[end], name = "h^3"),
    #scatter(x=hs,y=hs/hs[end]*err_EDPDL_field[end], name = "h")],
    Layout(
        yaxis_type="log", 
        xaxis_type="log",
        title = "Errors of the field for the exterior problem formualtions")
)