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
    ine::Nullable{HRepresentation{N, Float64}}
    ines::Nullable{MixedMatHRep{N, Float64}}
    ext::Nullable{VRepresentation{N, Float64}}
    exts::Nullable{MixedMatVRep{N, Float64}}
    noredundantinequality::Bool
    noredundantgenerator::Bool
    solver
    area::Nullable{Float64}
    volume::Nullable{Float64}

    function QHullPolyhedron{N}(ine, ext, nri::Bool, nrg::Bool, solver) where N
        new(ine, nothing, ext, nothing, nri, nrg, solver, nothing, nothing)
    end
    function QHullPolyhedron{N}(ine::HRepresentation{N, Float64}, solver) where N
        new(ine, nothing, nothing, nothing, false, false, solver, nothing, nothing)
    end
    function QHullPolyhedron{N}(ext::VRepresentation{N, Float64}, solver) where N
        new(nothing, nothing, ext, nothing, false, false, solver, nothing, nothing)
    end
end
Polyhedra.library(p::QHullPolyhedron) = QHullLibrary(p.solver)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, d::FullDim{1}, ::Type{Float64}) = default_type(d, Float64)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, d::FullDim, ::Type{T}) where T = default_type(d, T)
Polyhedra.similar_type(::Type{<:QHullPolyhedron}, ::FullDim{N}, ::Type{Float64}) where N = QHullPolyhedron{N}

function Polyhedra.similar(p::Tuple{<:QHullPolyhedron}, d::FullDim{1}, T::Type{Float64}, it::Polyhedra.It{1, Float64}...)
    Polyhedra.default_similar(p, d, T, it...)
end
function Polyhedra.similar(p::Tuple{<:QHullPolyhedron}, ::FullDim{N}, ::Type{Float64}, it::Polyhedra.It{N, Float64}...) where {N}
    QHullPolyhedron{N}(it..., p[1].solver)
end

function Polyhedra.arraytype(p::QHullPolyhedron)
    if isnull(p.ine) && !isnull(p.ines)
        p.ine = get(p.ines)
    end
    if isnull(p.ext) && !isnull(p.exts)
        p.ext = get(p.exts)
    end
    if isnull(p.ine)
        Polyhedra.arraytype(get(p.ext))
    elseif isnull(p.ext)
        Polyhedra.arraytype(get(p.ine))
    else
        @assert Polyhedra.arraytype(get(p.ine)) == Polyhedra.arraytype(get(p.ext))
        Polyhedra.arraytype(get(p.ine))
    end
end


# Helpers
epsz = 1e-8

function qhull(p::QHullPolyhedron{N}, rep=:Auto) where N
    if rep == :V || (rep == :Auto && (!isnull(p.ext) || !isnull(p.exts)))
        p.ext, ine, p.area, p.volume = qhull(getexts(p), p.solver)
        p.exts = nothing
        p.noredundantgenerator = true
        if isnull(p.ine)
            # Otherwise, it is not interesting as it may have redundancy
            p.ine = ine
            p.ines = nothing
        end
    else
        @assert rep == :H || rep == :Auto
        p.ine, ext, p.area, p.volume = qhull(getines(p), p.solver)
        p.ines = nothing
        p.noredundantinequality = true
        if isnull(p.ext)
            p.ext = ext
            p.exts = nothing
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
        if !isnull(p.ines)
            p.ine = p.ines
        else
            qhull(p)
        end
    end
    get(p.ine)
end
function getines(p::QHullPolyhedron)
    if isnull(p.ines)
        p.ines = MixedMatHRep(getine(p))
    end
    get(p.ines)
end
function getext(p::QHullPolyhedron)
    if isnull(p.ext)
        if !isnull(p.exts)
            p.ext = p.exts
        else
            qhull(p)
        end
    end
    get(p.ext)
end
function getexts(p::QHullPolyhedron)
    if isnull(p.exts)
        p.exts = MixedMatVRep(getext(p))
    end
    get(p.exts)
end

function clearfield!(p::QHullPolyhedron)
    p.ine = nothing
    p.ines = nothing
    p.ext = nothing
    p.exts = nothing
    p.noredundantinequality = false
    p.noredundantgenerator = false
end

# Implementation of Polyhedron's mandatory interface
function Polyhedra.polyhedron(repit::Representation{N}, lib::QHullLibrary) where {N}
    QHullPolyhedron{N}(repit, lib.solver)
end

QHullPolyhedron{N}(h::HRepresentation{N}, solver=nothing) where N = QHullPolyhedron{N}(MixedMatHRep{N,Float64}(h), solver)
QHullPolyhedron{N}(v::VRepresentation{N}, solver=nothing) where N = QHullPolyhedron{N}(MixedMatVRep{N,Float64}(v), solver)

QHullPolyhedron{N}(hps::Polyhedra.HyperPlaneIt{N}, hss::Polyhedra.HalfSpaceIt{N}, solver=nothing) where N = QHullPolyhedron{N}(MixedMatHRep{N, Float64}(hps, hss), solver)
QHullPolyhedron{N}(ps::Polyhedra.PointIt{N}, ls::Polyhedra.LineIt{N}, rs::Polyhedra.RayIt{N}, solver=nothing) where N = QHullPolyhedron{N}(MixedMatVRep{N, Float64}(ps, ls, rs), solver)

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
    p.ines = nothing
end
function Polyhedra.setvrep!(p::QHullPolyhedron{N}, v::VRepresentation{N}) where N
    p.ext = v
    p.exts = nothing
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
