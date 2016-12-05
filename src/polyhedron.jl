importall Base, Polyhedra
export QHullLib

type QHullLib <: PolyhedraLibrary
end

type QHullPolyhedron{N} <: Polyhedron{N, Float64}
    ine::Nullable{HRepresentation{N, Float64}}
    ines::Nullable{SimpleHRepresentation{N, Float64}}
    ext::Nullable{VRepresentation{N, Float64}}
    exts::Nullable{SimpleVRepresentation{N, Float64}}
    hlinearitydetected::Bool
    vlinearitydetected::Bool
    noredundantinequality::Bool
    noredundantgenerator::Bool

    function QHullPolyhedron(ine::HRepresentation{N, Float64}, ext::VRepresentation{N, Float64}, hld::Bool, vld::Bool, nri::Bool, nrg::Bool)
        new(ine, nothing, ext, nothing, hld, vld, nri, nrg)
    end
    function QHullPolyhedron(ine::HRepresentation{N, Float64})
        new(ine, nothing, nothing, nothing, false, false, false, false)
    end
    function QHullPolyhedron(ext::VRepresentation{N, Float64})
        new(nothing, nothing, ext, nothing, false, false, false, false)
    end
end

# ine may decompose fast but if ine is nothing I do not want to ask to compute it to see the type it is
# saying false normally do not give troubles
decomposedhfast{N}(::Type{QHullPolyhedron{N}}) = false
decomposedvfast{N}(::Type{QHullPolyhedron{N}}) = false
decomposedhfast{N}(p::QHullPolyhedron{N}) = decomposedhfast(QHullPolyhedron{N})
decomposedvfast{N}(p::QHullPolyhedron{N}) = decomposedvfast(QHullPolyhedron{N})

eltype{N}(::Type{QHullPolyhedron{N}}) = Float64
eltype(::QHullPolyhedron) = Float64

# Helpers
epsz = 1e-8

function qhull{N}(p::QHullPolyhedron{N}, rep=:Auto)
    if rep == :V || (rep == :Auto && (!isnull(p.exts) || !isnull(p.exts)))
        p.ext, ine = qhull(getexts(p))
        p.exts = nothing
        p.noredundantgenerator = true
        if isnull(p.ine)
            # Otherwise, it is not interesting as it may have redundancy
            p.ine = ine
            p.ines = nothing
        end
    else
        @assert rep == :H || rep == :Auto
        p.ine, ext = qhull(getines(p))
        p.ines = nothing
        p.noredundantinequality = true
        if isnull(p.ext)
            p.ext = ext
            p.exts = nothing
        end
    end
end

function qhull{N, T}(h::SimpleVRepresentation{N, T})
    if hasrays(h)
        error("Rays are not supported.")
    end
    V = h.V
    linset = h.Vlinset
    clinset = collect(linset)
    if !isempty(linset)
        V = [V; -(@view V[clinset, :])]
    end
    ch = chull(V)
    V = ch.points[ch.vertices, :]
    vnored = SimpleVRepresentation(V)
    A = ch.facets[:, 1:N]
    b = ch.facets[:, N+1]
    h = SimpleHRepresentation(A, b)
    vnored, h
end

function qhull{N, T<:Real}(h::SimpleHRepresentation{N, T})
    A = h.A
    b = h.b
    linset = h.linset
    if !isempty(linset)
        error("The origin should be in the interior of the polytope so equalities are not supported.")
    end
    m = size(A, 1)
    B = Matrix{T}(m, N)
    for i in 1:m
        if i in linset || b[i] < epsz
            error("The origin should be in the interior of the polytope but the $(i)th inequality is not safisfied at the origin.")
        end
        B[i,:] = (@view A[i,:]) / b[i]
    end
    ch = chull(B)
    C = ch.points[ch.vertices, :]
    hnored = SimpleHRepresentation(C, ones(size(C, 1)))
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
            R[nr, :] = @view Vlift[i, 1:N]
        else
            nv += 1
            V[nv, :] = (@view Vlift[i, 1:N]) / Vlift[i, N+1]
        end
    end
    v = SimpleVRepresentation(V, R)
    hnored, v
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
        p.ines = SimpleHRepresentation(getine(p))
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
        p.exts = SimpleVRepresentation(getext(p))
    end
    get(p.exts)
end

function clearfield!(p::QHullPolyhedron)
    p.ine = nothing
    p.ines = nothing
    p.ext = nothing
    p.exts = nothing
    hlinearitydetected = false
    vlinearitydetected = false
    noredundantinequality = false
    noredundantgenerator = false
end
function updateine!{N}(p::QHullPolyhedron{N}, ine::HRepresentation{N, Float64})
    clearfield!(p)
    p.ine = ine
end
function updateext!{N}(p::QHullPolyhedron{N}, ext::VRepresentation{N, Float64})
    clearfield!(p)
    p.ext = ext
end


# Implementation of Polyhedron's mandatory interface
polyhedron{N}(repit::Union{Representation{N},HRepIterator{N},VRepIterator{N}}, ::QHullLib) = QHullPolyhedron{N}(repit)

getlibraryfor{T<:AbstractFloat}(::QHullPolyhedron, ::Type{T}) = QHullLib()

(::Type{QHullPolyhedron{N}}){N, T}(it::HRepIterator{N,T}) = QHullPolyhedron{N}(SimpleHRepresentation{N,Float64}(it))
(::Type{QHullPolyhedron{N}}){N, T}(it::VRepIterator{N,T}) = QHullPolyhedron{N}(SimpleVRepresentation{N,Float64}(it))

function (::Type{QHullPolyhedron{N}}){N}(; eqs=nothing, ineqs=nothing, points=nothing, rays=nothing)
    noth = eqs === nothing && ineqs === nothing
    notv = points === nothing && rays === nothing
    if noth && notv
        error("QHullPolyhedron should have at least one iterator to be built")
    end
    if !noth && !notv
        error("QHullPolyhedron constructed with a combination of eqs/ineqs with points/rays")
    end
    if notv
        QHullPolyhedron{N}(SimpleHRepresentation{N, Float64}(eqs=eqs, ineqs=ineqs))
    else
        QHullPolyhedron{N}(SimpleVRepresentation{N, Float64}(points=points, rays=rays))
    end
end

function Base.copy{N}(p::QHullPolyhedron{N})
    ine = nothing
    if !isnull(p.ine)
        ine = copy(get(p.ine))
    end
    ext = nothing
    if !isnull(p.ext)
        ext = copy(get(p.ext))
    end
    QHullPolyhedron{N}(ine, ext, p.hlinearitydetected, p.vlinearitydetected, p.noredundantinequality, p.noredundantgenerator)
end
function Base.push!{N}(p::QHullPolyhedron{N}, ine::HRepresentation{N})
    updateine!(p, intersect(getine(p), changeeltype(ine, Float64)))
end
function Base.push!{N}(p::QHullPolyhedron{N}, ext::VRepresentation{N})
    updateext!(p, getext(p) + changeeltype(ext, Float64))
end
function hrepiscomputed(p::QHullPolyhedron)
    !isnull(p.ine)
end
function hrep(p::QHullPolyhedron)
    copy(getine(p))
end
function vrepiscomputed(p::QHullPolyhedron)
    !isnull(p.ext)
end
function vrep(p::QHullPolyhedron)
    copy(getext(p))
end
function detecthlinearities!(p::QHullPolyhedron)
    warn("detecthlinearities! not supported yet")
end
function detectvlinearities!(p::QHullPolyhedron)
    warn("detectvlinearities! not supported yet")
end
function removehredundancy!(p::QHullPolyhedron)
    qhull(p, :H)
end
function removevredundancy!(p::QHullPolyhedron)
    qhull(p, :V)
end

for f in [:hashreps, :nhreps, :starthrep, :hasineqs, :nineqs, :startineq, :haseqs, :neqs, :starteq]
    @eval $f(p::QHullPolyhedron) = $f(getine(p))
end
for f in [:donehrep, :nexthrep, :doneineq, :nextineq, :doneeq, :nexteq]
    @eval $f(p::QHullPolyhedron, state) = $f(getine(p), state)
end

for f in [:hasvreps, :nvreps, :startvrep, :haspoints, :npoints, :startpoint, :hasrays, :nrays, :startray]
    @eval $f(p::QHullPolyhedron) = $f(getext(p))
end
for f in [:donevrep, :nextvrep, :donepoint, :nextpoint, :doneray, :nextray]
    @eval $f(p::QHullPolyhedron, state) = $f(getext(p), state)
end
