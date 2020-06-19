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
    hrep::Union{Nothing, Polyhedra.MixedMatHRep{Float64, Matrix{Float64}}}
    vrep::Union{Nothing, Polyhedra.MixedMatVRep{Float64, Matrix{Float64}}}
    nohredundancy::Bool
    novredundancy::Bool
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

Polyhedra.FullDim(p::Polyhedron) = Polyhedra.FullDim_rep(p.hrep, p.vrep)
Polyhedra.library(p::Polyhedron) = Library(p.solver)
Polyhedra.similar_type(::Type{<:Polyhedron}, d::Polyhedra.FullDim, T::Type) = Polyhedra.default_type(d, T)
function Polyhedra.similar_type(::Type{<:Polyhedron}, d::Polyhedra.FullDim, ::Type{Float64})
    if _isone(d)
        return Polyhedra.default_type(d, Float64)
    else
        Polyhedron
    end
end

Polyhedra.default_solver(p::Polyhedron; T=nothing) = p.solver
Polyhedra.supportssolver(::Type{<:Polyhedron}) = true

Polyhedra.hvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.hvectortype(Polyhedra.MixedMatHRep{Float64, Matrix{Float64}})
Polyhedra.vvectortype(::Union{Polyhedron, Type{Polyhedron}}) = Polyhedra.vvectortype(Polyhedra.MixedMatVRep{Float64, Matrix{Float64}})

# Helpers
epsz = 1e-8

function qhull(p::Polyhedron, rep=:Auto)
    if rep == :V || (rep == :Auto && (p.vrep !== nothing))
        p.vrep, ine, p.area, p.volume = qhull(Polyhedra.vrep(p), p.solver)
        if p.hrep === nothing
            # Otherwise, it is not interesting as it may have redundancy
            p.hrep = ine
        end
    else
        @assert rep == :H || rep == :Auto
        p.hrep, ext = qhull(Polyhedra.hrep(p), p.solver)
        if p.vrep === nothing
            p.vrep = ext
        end
    end
    p.nohredundancy = true
    p.novredundancy = true
    return
end

function qhull(h::Polyhedra.MixedMatVRep{T}, solver) where T
    if Polyhedra.hasrays(h)
        error("Rays are not supported.")
    end
    V = h.V
    ch = chull(V)
    V = ch.points[ch.vertices, :]
    vnored = Polyhedra.vrep(V)
    N = Polyhedra.fulldim(h)
    A = ch.facets[:, 1:N]
    b = -ch.facets[:, N+1]
    h = Polyhedra.hrep(A, b)
    # The same facet may be listed several times in case the facet is not
    # simplicial, e.g. a quad facet is listed twice as it is made of two
    # triangles. After `removeduplicates`, there should be no redundancy.
    hnored = Polyhedra.removeduplicates(h, Polyhedra.OppositeMockOptimizer)
    vnored, hnored, ch.area, ch.volume
end

function qhull(h::Polyhedra.MixedMatHRep{T}, solver) where {T<:Real}
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
    # See comment in `qhull(::MixedMatVRep)`.
    vnored = Polyhedra.removeduplicates(v)
    # `ch.area`, `ch.volume` give area and volume of the polar of the shifted
    # polytope, not the area and volume of the polytope with reps `h` and `v`.
    hnored, vnored
end

function clearfield!(p::Polyhedron)
    p.hrep = nothing
    p.vrep = nothing
    p.nohredundancy = false
    p.novredundancy = false
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
    if p.hrep !== nothing
        ine = copy(p.hrep)
    end
    ext = nothing
    if p.vrep !== nothing
        ext = copy(p.vrep)
    end
    Polyhedron(ine, ext, p.nohredundancy, p.novredundancy, p.solver)
end
function Polyhedra.removehredundancy!(p::Polyhedron)
    if !p.nohredundancy
        qhull(p, :H)
    end
end
function Polyhedra.removevredundancy!(p::Polyhedron)
    if !p.novredundancy
        qhull(p, :V)
    end
end

function Polyhedra.volume(p::Polyhedron)
    if p.volume === nothing
        qhull(p, :V)
    end
    return p.volume
end

Polyhedra.hrepiscomputed(p::Polyhedron) = p.hrep !== nothing
Polyhedra.computehrep!(p::Polyhedron) = qhull(p, :V)
function Polyhedra.hrep(p::Polyhedron)
    if !Polyhedra.hrepiscomputed(p)
        Polyhedra.computehrep!(p)
    end
    return p.hrep
end
Polyhedra.vrepiscomputed(p::Polyhedron) = p.vrep !== nothing
Polyhedra.computevrep!(p::Polyhedron) = qhull(p, :H)
function Polyhedra.vrep(p::Polyhedron)
    if !Polyhedra.vrepiscomputed(p)
        Polyhedra.computevrep!(p)
    end
    return p.vrep
end

function Polyhedra.sethrep!(p::Polyhedron, h::Polyhedra.HRepresentation)
    p.hrep = h
end
function Polyhedra.setvrep!(p::Polyhedron, v::Polyhedra.VRepresentation)
    p.vrep = v
end
function Polyhedra.resethrep!(p::Polyhedron, h::Polyhedra.HRepresentation)
    clearfield!(p)
    p.hrep = h
end
function Polyhedra.resetvrep!(p::Polyhedron, v::Polyhedra.VRepresentation)
    clearfield!(p)
    p.vrep = v
end
