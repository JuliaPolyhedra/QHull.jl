import StaticArrays
import Polyhedra

struct Library <: Polyhedra.Library
    solver
end
Library() = Library(nothing)
_isone(d::Int) = isone(d)
_isone(d::StaticArrays.Size{(1,)}) = true
_isone(d::StaticArrays.Size) where T = false
Polyhedra.similar_library(l::Library, d::Polyhedra.FullDim, T::Type) = Polyhedra.default_library(d, T)
function Polyhedra.similar_library(l::Library, d::Polyhedra.FullDim, ::Type{Float64})
    if _isone(d)
        return Polyhedra.default_library(d, Float64)
    else
        return Library(l.solver)
    end
end

mutable struct Polyhedron <: Polyhedra.Polyhedron{Float64}
    ine::Union{Nothing, Polyhedra.MixedMatHRep{Float64, Matrix{Float64}}}
    ext::Union{Nothing, Polyhedra.MixedMatVRep{Float64, Matrix{Float64}}}
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver
    area::Union{Nothing, Float64}
    volume::Union{Nothing, Float64}

    function Polyhedron(ine, ext, nri::Bool, nrg::Bool, solver)
        new(ine, ext, nri, nrg, solver, nothing, nothing)
    end
    function Polyhedron(ine::Polyhedra.HRepresentation{Float64}, solver)
        new(ine, nothing, false, false, solver, nothing, nothing)
    end
    function Polyhedron(ext::Polyhedra.VRepresentation{Float64}, solver)
        new(nothing, ext, false, false, solver, nothing, nothing)
    end
end

Polyhedra.FullDim(p::Polyhedron) = Polyhedra.FullDim_rep(p.ine, p.ext)
Polyhedra.library(p::Polyhedron) = Library(p.solver)
Polyhedra.similar_type(::Type{<:Polyhedron}, d::Polyhedra.FullDim, T::Type) = Polyhedra.default_type(d, T)
function Polyhedra.similar_type(::Type{<:Polyhedron}, d::Polyhedra.FullDim, ::Type{Float64})
    if _isone(d)
        return Polyhedra.default_type(d, Float64)
    else
        Polyhedron
    end
end

Polyhedra.default_solver(p::Polyhedron) = p.solver
Polyhedra.supportssolver(::Type{<:Polyhedron}) = true

Polyhedra.hvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.hvectortype(Polyhedra.MixedMatHRep{Float64, Matrix{Float64}})
Polyhedra.vvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.vvectortype(Polyhedra.MixedMatVRep{Float64, Matrix{Float64}})

# Helpers
epsz = 1e-8

function qhull(p::Polyhedron, rep=:Auto)
    if rep == :V || (rep == :Auto && (p.ext !== nothing))
        p.ext, ine, p.area, p.volume = qhull(getext(p), p.solver)
        p.noredundantgenerator = true
        if p.ine === nothing
            # Otherwise, it is not interesting as it may have redundancy
            p.ine = ine
        end
    else
        @assert rep == :H || rep == :Auto
        p.ine, ext, p.area, p.volume = qhull(getine(p), p.solver)
        p.noredundantinequality = true
        if p.ext === nothing
            p.ext = ext
        end
    end
end

function qhull(h::Polyhedra.MixedMatVRep{T}, solver=nothing) where T
    if Polyhedra.hasrays(h)
        error("Rays are not supported.")
    end
    V = h.V
    ch = chull(V)
    V = ch.points[ch.vertices, :]
    vnored = Polyhedra.vrep(V)
    N = Polyhedra.fulldim(h)
    A = ch.facets[:, 1:N]
    b = ch.facets[:, N+1]
    h = Polyhedra.hrep(A, b)
    vnored, h, ch.area, ch.volume
end

function qhull(h::Polyhedra.MixedMatHRep{T}, solver=nothing) where {T<:Real}
    linset = h.linset
    if !isempty(linset)
        error("Equalities are not supported.")
    end
    N = Polyhedra.fulldim(h)
    containorigin = Polyhedra.ininterior(zeros(N), h)
    if !containorigin
        chebycenter, chebyradius = Polyhedra.chebyshevcenter(h, solver)
        h = Polyhedra.translate(h, -chebycenter)
    end

    A = h.A
    b = h.b
    m = size(A, 1)
    B = Matrix{T}(undef, m, N)
    for i in 1:m
        @assert !(i in linset)
        if b[i] < epsz
            error("The origin should be in the interior of the polytope but the $(i)th inequality is not safisfied at the origin.")
        end
        B[i,:] = (@view A[i,:]) / b[i]
    end
    ch = chull(B)
    C = ch.points[ch.vertices, :]
    hnored = Polyhedra.hrep(C, ones(size(C, 1)))
    Vlift = ch.facets
    nvreps = size(Vlift, 1)
    irays = BitSet()
    ipoints = BitSet()
    for i in 1:nvreps
        if Vlift[i, N+1] > epsz
            error("The origin should be in the interior of the convex hull")
        end
        if Vlift[i, N+1] > -epsz
            push!(irays, i)
        else
            push!(ipoints, i)
        end
    end
    rays = collect(irays)
    points = collect(ipoints)
    R = Matrix{T}(undef, length(rays), N)
    V = Matrix{T}(undef, length(points), N)
    nr = nv = 0
    for i in 1:nvreps
        if i in irays
            nr += 1
            R[nr, :] = -@view Vlift[i, 1:N]
        else
            nv += 1
            V[nv, :] = -(@view Vlift[i, 1:N]) / Vlift[i, N+1]
        end
    end
    v = Polyhedra.vrep(V, R)
    if !containorigin
        hnored = Polyhedra.translate(hnored, chebycenter)
    end
    if !containorigin
        v = Polyhedra.translate(v, chebycenter)
    end
    hnored, v, ch.area, ch.volume
end

function getine(p::Polyhedron)
    if p.ine === nothing
        qhull(p)
    end
    p.ine
end
function getext(p::Polyhedron)
    if p.ext === nothing
        qhull(p)
    end
    p.ext
end

function clearfield!(p::Polyhedron)
    p.ine = nothing
    p.ext = nothing
    p.noredundantinequality = false
    p.noredundantgenerator = false
end

# Implementation of Polyhedron's mandatory interface
function Polyhedra.polyhedron(repit::Polyhedra.Representation, lib::Library)
    Polyhedron(repit, lib.solver)
end

Polyhedron(h::Polyhedra.HRepresentation, solver) = Polyhedron(convert(Polyhedra.MixedMatHRep{Float64, Matrix{Float64}}, h), solver)
Polyhedron(v::Polyhedra.VRepresentation, solver) = Polyhedron(convert(Polyhedra.MixedMatVRep{Float64, Matrix{Float64}}, v), solver)

Polyhedron(d::Polyhedra.FullDim, hps::Polyhedra.HyperPlaneIt, hss::Polyhedra.HalfSpaceIt; solver=nothing) = Polyhedron(Polyhedra.MixedMatHRep{Float64, Matrix{Float64}}(d, hps, hss), solver)
Polyhedron(d::Polyhedra.FullDim, ps::Polyhedra.PointIt, ls::Polyhedra.LineIt, rs::Polyhedra.RayIt; solver=nothing) = Polyhedron(Polyhedra.MixedMatVRep{Float64, Matrix{Float64}}(d, ps, ls, rs), solver)

function Base.copy(p::Polyhedron)
    ine = nothing
    if p.ine !== nothing
        ine = copy(p.ine)
    end
    ext = nothing
    if p.ext !== nothing
        ext = copy(p.ext)
    end
    Polyhedron(ine, ext, p.noredundantinequality, p.noredundantgenerator, p.solver)
end
function Polyhedra.hrepiscomputed(p::Polyhedron)
    p.ine !== nothing
end
function Polyhedra.hrep(p::Polyhedron)
    getine(p)
end
function Polyhedra.vrepiscomputed(p::Polyhedron)
    p.ext !== nothing
end
function Polyhedra.vrep(p::Polyhedron)
    getext(p)
end
function Polyhedra.sethrep!(p::Polyhedron, h::Polyhedra.HRepresentation)
    p.ine = h
end
function Polyhedra.setvrep!(p::Polyhedron, v::Polyhedra.VRepresentation)
    p.ext = v
end
function Polyhedra.resethrep!(p::Polyhedron, h::Polyhedra.HRepresentation)
    clearfield!(p)
    p.ine = h
end
function Polyhedra.resetvrep!(p::Polyhedron, v::Polyhedra.VRepresentation)
    clearfield!(p)
    p.ext = v
end
function Polyhedra.removehredundancy!(p::Polyhedron)
    qhull(p, :H)
end
function Polyhedra.removevredundancy!(p::Polyhedron)
    qhull(p, :V)
end

function Polyhedra.volume(p::Polyhedron)
    if p.volume === nothing
        qhull(p)
    end
    p.volume
end
