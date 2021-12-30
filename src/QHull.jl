## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module QHull

export Chull, chull

using PyCall
const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end

mutable struct Chull{T<:Real}
    points::Matrix{T}
    vertices::Vector{Int}
    simplices::Vector{Vector{Int}}
    facets::Matrix{T}
    area::Float64
    volume::Float64
end

## helper for base-0 / base-1 difference
incone(x) = for i in 1:length(x)
    x[i] += 1
end

import QHull_jll

function chull(x::Matrix{T}) where T<:Real
    q = ccall((:qh_alloc_qh, Qhull_jll.libqhull_r), Ptr{Cvoid}, (Ptr{Cvoid},), Base.stderr)
    ccall((:qh_init_B, Qhull_jll.libqhull_r), Cvoid, (Ptr{Cvoid}, Ptr{Cdouble}, Cint, Cint, Cuint), q, Matrix(x'), size(x, 1), size(x, 2), 0)
    py = spatial.ConvexHull(x)
    points = convert(Matrix{T}, py."points")
    vertices = convert(Vector{Int}, py."vertices")
    incone(vertices)
    simplices = convert(Vector{Vector{Int}}, py."simplices")
    for simplex in simplices
        incone(simplex)
    end
    facets = convert(Matrix{T}, py."equations")
    area = convert(Float64, py."area")
    volume = convert(Float64, py."volume")
    Chull(points, vertices, simplices, facets, area, volume)
end

include("polyhedron.jl")

function Base.show(io::IO, ::MIME"text/plain", ch::Chull)
    println(io, string("Convex Hull of ", size(ch.points, 1), " points in ", size(ch.points, 2), " dimensions"))
    println(io, "Hull segment vertex indices:")
    println(io, ch.vertices)
    println(io, "Points on convex hull in original order:\n")
    println(io, ch.points[sort(ch.vertices[:, 1]), :])
end

end
