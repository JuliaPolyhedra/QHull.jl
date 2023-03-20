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
    vertices::Vector{Int32}
    simplices::Matrix{Int32}
    facets::Matrix{T}
    area::Float64
    volume::Float64
end

function chull(x::Matrix{T}) where T<:Real
    py = spatial.ConvexHull(x)
    res = Chull(py.points, py.vertices, py.simplices, py.equations, py.area, py.volume)
    res.vertices .+= 1
    res.simplices .+= 1
    res
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
