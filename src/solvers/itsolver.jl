import IterativeSolvers



struct GMRESSolver{L,T} <: LinearMap{T}
    linear_operator::L
    maxiter::Int
    restart::Int
    tol::T
    verbose::Bool
end


Base.axes(A::GMRESSolver) = reverse(axes(A.linear_operator))


function GMRESSolver(op; maxiter=0, restart=0, tol=sqrt(eps(real(eltype(op)))), verbose=true)

    m, n = size(op)
    @assert m == n

    maxiter == 0 && (maxiter = div(n, 5))
    restart == 0 && (restart = n)

    GMRESSolver(op, maxiter, restart, tol, verbose)
end


operator(solver::GMRESSolver) = solver.linear_operator


function solve(solver::GMRESSolver, b; abstol=zero(real(eltype(b))), reltol=solver.tol)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Array{T}, axes(solver)[2])
    fill!(x,0)
    x, ch = solve!(x, solver, b; abstol, reltol)
end


function solve!(x, solver::GMRESSolver, b; abstol=zero(real(eltype(b))), reltol=solver.tol)
    op = operator(solver)
    x, ch = IterativeSolvers.gmres!(x, op, b, log=true,  maxiter=solver.maxiter,
        restart=solver.restart, reltol=reltol, abstol=abstol, verbose=solver.verbose)
    return x, ch
end


function Base.:*(A::GMRESSolver, b::AbstractVector)

    T = promote_type(eltype(A), eltype(b))
    y = PseudoBlockVector{T}(undef, (axes(A,2),))

    mul!(y, A, b)
end

Base.size(solver::GMRESSolver) = reverse(size(solver.linear_operator))

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::GMRESSolver, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    solver.verbose && println("Number of iterations: ", ch.iters)
    ch.isconverged || error("Iterative solver did not converge.")
    return y
end

LinearAlgebra.adjoint(A::GMRESSolver) = GMRESSolver(adjoint(A.linear_operator), A.maxiter, A.restart, A.tol, A.verbose)
LinearAlgebra.transpose(A::GMRESSolver) = GMRESSolver(transpose(A.linear_operator), A.maxiter, A.restart, A.tol, A.verbose)



function gmres_ch(eq::DiscreteEquation; maxiter=0, restart=0, tol=0)

    lhs = eq.equation.lhs
    rhs = eq.equation.rhs

    X = _spacedict_to_directproductspace(eq.test_space_dict)
    Y = _spacedict_to_directproductspace(eq.trial_space_dict)

    b = assemble(rhs, X)
    Z = assemble(lhs, X, Y)

    if tol == 0
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart)
    else
        invZ = GMRESSolver(Z, maxiter=maxiter, restart=restart, tol=tol)
    end
    x, ch = solve(invZ, b)
    # x = invZ * b

    ax = nestedrange(Y, 1, numfunctions)
    return PseudoBlockVector(x, (ax,)), ch
end

gmres(eq::DiscreteEquation; maxiter=0, restart=0, tol=0) = gmres_ch(eq; maxiter, restart, tol)[1]



struct CGSolver{L,M,T} <: LinearMap{T}
    A::L
    Pl::M
    abstol::T
    reltol::T
    maxiter::Int
    verbose::Bool
end


Base.size(solver::CGSolver) = reverse(size(solver.A))
Base.axes(solver::CGSolver) = reverse(axes(solver.A))
operator(solver::CGSolver) = solver.A


cg(A;
    Pl = IterativeSolvers.Identity(),
    abstol::Real = zero(real(eltype(A))),
    reltol::Real = sqrt(eps(real(eltype(A)))),
    maxiter::Int = size(A,2),
    verbose::Bool = false) = CGSolver(A, Pl, abstol, reltol, maxiter, verbose)


function solve(solver::CGSolver, b)
    T = promote_type(eltype(solver), eltype(b))
    x = similar(Vector{T}, size(solver)[2])
    fill!(x,0)
    x, ch = solve!(x, solver, b)
    z = similar(Array{T}, axes(solver)[2])
    copyto!(z,x)
    return z, ch
end


function solve!(x, solver::CGSolver, b)
    op = operator(solver)
    x, ch = IterativeSolvers.cg!(x, op, b;
        Pl = solver.Pl,
        abstol = solver.abstol,
        reltol = solver.reltol,
        maxiter = solver.maxiter,
        verbose = solver.verbose,
        log = true)
    return x, ch
end

function LinearAlgebra.mul!(y::AbstractVecOrMat, solver::CGSolver, x::AbstractVector)
    fill!(y,0)
    y, ch = solve!(y, solver, x)
    return y
end


struct Preconditioner{L,T} <: LinearMap{T}
    A::L
end

Preconditioner(A::L) where {L} = Preconditioner{L,eltype(A)}(A)

function LinearAlgebra.ldiv!(y, P::Preconditioner, x)
    mul!(y, P.A, x)
end