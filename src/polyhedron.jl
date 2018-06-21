using Polyhedra
export QHullLibrary
using MathProgBase

struct QHullLibrary <: PolyhedraLibrary
    solver
end
QHullLibrary() = QHullLibrary(nothing)
Polyhedra.similar_library(l::QHullLibrary, d::FullDim{1}, ::Type{T}) where T = default_library(d, T)
Polyhedra.similar_library(l::QHullLibrary, d::FullDim, ::Type{T}) where T = default_library(d, T)
Polyhedra.similar_library(l::QHullLibrary, ::FullDim, ::Type{<:AbstractFloat}) = QHullLibrary(l.solver)

mutable struct QHullPolyhedron{N} <: Polyhedron{N, Float64}
    ine::Nullable{MixedMatHRep{N, Float64, Matrix{Float64}}}
    ext::Nullable{MixedMatVRep{N, Float64, Matrix{Float64}}}
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver
    area::Nullable{Float64}
    volume::Nullable{Float64}

    function QHullPolyhedron{N}(ine, ext, nri::Bool, nrg::Bool, solver) where N
        new(ine, ext, nri, nrg, solver, nothing, nothing)
    end
    function QHullPolyhedron{N}(ine::HRepresentation{N, Float64}, solver) where N
        new(ine, nothing, false, false, solver, nothing, nothing)
    end
    function QHullPolyhedron{N}(ext::VRepresentation{N, Float64}, solver) where N
        new(nothing, ext, false, false, solver, nothing, nothing)
    end
end
Polyhedra.library(p::QHullPolyhedron) = QHullLibrary(p.solver)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, d::FullDim{1}, ::Type{Float64}) = default_type(d, Float64)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, d::FullDim, ::Type{T}) where T = default_type(d, T)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, ::FullDim{N}, ::Type{Float64}) where N = QHullPolyhedron{N}

Polyhedra.default_solver(p::QHullPolyhedron) = p.solver
Polyhedra.supportssolver(::Type{<:QHullPolyhedron}) = true

Polyhedra.hvectortype(::Union{QHullPolyhedron{N}, Type{QHullPolyhedron{N}}}) where N = Polyhedra.hvectortype(MixedMatHRep{N, Float64, Matrix{Float64}})
Polyhedra.vvectortype(::Union{QHullPolyhedron{N}, Type{QHullPolyhedron{N}}}) where N = Polyhedra.vvectortype(MixedMatVRep{N, Float64, Matrix{Float64}})

# Helpers
epsz = 1e-8

function qhull(p::QHullPolyhedron{N}, rep=:Auto) where N
    if rep == :V || (rep == :Auto && (!isnull(p.ext)))
        p.ext, ine, p.area, p.volume = qhull(getext(p), p.solver)
        p.noredundantgenerator = true
        if isnull(p.ine)
            # Otherwise, it is not interesting as it may have redundancy
            p.ine = ine
        end
    else
        @assert rep == :H || rep == :Auto
        p.ine, ext, p.area, p.volume = qhull(getine(p), p.solver)
        p.noredundantinequality = true
        if isnull(p.ext)
            p.ext = ext
        end
    end
end

function qhull(h::MixedMatVRep{N, T}, solver=nothing) where {N, T}
    if hasrays(h)
        error("Rays are not supported.")
    end
    V = h.V
    ch = chull(V)
    V = ch.points[ch.vertices, :]
    vnored = vrep(V)
    A = ch.facets[:, 1:N]
    b = ch.facets[:, N+1]
    h = hrep(A, b)
    vnored, h, ch.area, ch.volume
end

function qhull(h::MixedMatHRep{N, T}, solver=nothing) where {N, T<:Real}
    linset = h.linset
    if !isempty(linset)
        error("Equalities are not supported.")
    end
    containorigin = ininterior(zeros(N), h)
    if !containorigin
        chebycenter, chebyradius = chebyshevcenter(h, solver)
        h = translate(h, -chebycenter)
    end

    A = h.A
    b = h.b
    m = size(A, 1)
    B = Matrix{T}(m, N)
    for i in 1:m
        @assert !(i in linset)
        if b[i] < epsz
            error("The origin should be in the interior of the polytope but the $(i)th inequality is not safisfied at the origin.")
        end
        B[i,:] = (@view A[i,:]) / b[i]
    end
    ch = chull(B)
    C = ch.points[ch.vertices, :]
    hnored = hrep(C, ones(size(C, 1)))
    Vlift = ch.facets
    nvreps = size(Vlift, 1)
    irays = IntSet()
    ipoints = IntSet()
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
    R = Matrix{T}(length(rays), N)
    V = Matrix{T}(length(points), N)
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
    v = vrep(V, R)
    if !containorigin
        hnored = translate(hnored, chebycenter)
    end
    if !containorigin
        v = translate(v, chebycenter)
    end
    hnored, v, ch.area, ch.volume
end

function getine(p::QHullPolyhedron)
    if isnull(p.ine)
        qhull(p)
    end
    get(p.ine)
end
function getext(p::QHullPolyhedron)
    if isnull(p.ext)
        qhull(p)
    end
    get(p.ext)
end

function clearfield!(p::QHullPolyhedron)
    p.ine = nothing
    p.ext = nothing
    p.noredundantinequality = false
    p.noredundantgenerator = false
end

# Implementation of Polyhedron's mandatory interface
function Polyhedra.polyhedron(repit::Representation{N}, lib::QHullLibrary) where {N}
    QHullPolyhedron{N}(repit, lib.solver)
end

QHullPolyhedron{N}(h::HRepresentation{N}, solver) where N = QHullPolyhedron{N}(MixedMatHRep{N,Float64, Matrix{Float64}}(h), solver)
QHullPolyhedron{N}(v::VRepresentation{N}, solver) where N = QHullPolyhedron{N}(MixedMatVRep{N,Float64, Matrix{Float64}}(v), solver)

QHullPolyhedron{N}(hps::Polyhedra.HyperPlaneIt{N}, hss::Polyhedra.HalfSpaceIt{N}; solver=nothing) where N = QHullPolyhedron{N}(MixedMatHRep{N, Float64, Matrix{Float64}}(hps, hss), solver)
QHullPolyhedron{N}(ps::Polyhedra.PointIt{N}, ls::Polyhedra.LineIt{N}, rs::Polyhedra.RayIt{N}; solver=nothing) where N = QHullPolyhedron{N}(MixedMatVRep{N, Float64, Matrix{Float64}}(ps, ls, rs), solver)

function Base.copy(p::QHullPolyhedron{N}) where N
    ine = nothing
    if !isnull(p.ine)
        ine = copy(get(p.ine))
    end
    ext = nothing
    if !isnull(p.ext)
        ext = copy(get(p.ext))
    end
    QHullPolyhedron{N}(ine, ext, p.noredundantinequality, p.noredundantgenerator, p.solver)
end
function Polyhedra.hrepiscomputed(p::QHullPolyhedron)
    !isnull(p.ine)
end
function Polyhedra.hrep(p::QHullPolyhedron)
    getine(p)
end
function Polyhedra.vrepiscomputed(p::QHullPolyhedron)
    !isnull(p.ext)
end
function Polyhedra.vrep(p::QHullPolyhedron)
    getext(p)
end
function Polyhedra.sethrep!(p::QHullPolyhedron{N}, h::HRepresentation{N}) where N
    p.ine = h
end
function Polyhedra.setvrep!(p::QHullPolyhedron{N}, v::VRepresentation{N}) where N
    p.ext = v
end
function resethrep!(p::QHullPolyhedron{N}, h::HRepresentation{N}) where N
    clearfield!(p)
    p.ine = h
end
function resetvrep!(p::QHullPolyhedron{N}, v::VRepresentation{N}) where N
    clearfield!(p)
    p.ext = v
end
function Polyhedra.removehredundancy!(p::QHullPolyhedron)
    qhull(p, :H)
end
function Polyhedra.removevredundancy!(p::QHullPolyhedron)
    qhull(p, :V)
end

function Polyhedra.volume(p::QHullPolyhedron)
    if isnull(p.volume)
        qhull(p)
    end
    get(p.volume)
end
