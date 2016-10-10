## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

__precompile__(true)
module QHull

export Chull, chull, display, show

using PyCall
const spatial = PyNULL()

function __init__()
    copy!(spatial, pyimport_conda("scipy.spatial", "scipy"))
end

type Chull{T<:Real}
    points::Array{T}
    vertices::Vector{Int}
    simplices::Vector{Vector{Int}}
end

## helper for base-0 / base-1 difference
incone(x) = for i in 1:length(x)
    x[i] += 1
end

function chull{T<:Real}(x::Array{T})
    r, c = size(x)
    py = spatial[:ConvexHull](x)
    points = convert(Array{T}, py["points"])
    vertices = convert(Array{Int}, py["vertices"])
    incone(vertices)
    simplices = convert(Vector{Vector{Int}}, py["simplices"])
    for simplex in simplices
        incone(simplex)
    end
    Chull(points, vertices, simplices)
end

function Base.show(io::IO, ::MIME"text/plain", ch::Chull)
    println(io, string("Convex Hull of ", size(ch.points, 1), " points in ", size(ch.points, 2), " dimensions"))
    println(io, "Hull segment vertex indices:")
    println(io, ch.vertices)
    println(io, "Points on convex hull in original order:\n")
    println(io, ch.points[sort(ch.vertices[:, 1]), :])
end

using RecipesBase
@recipe function f{T<:Chull}(val::T)
    size(val.points, 2) > 2 && warning("Only the two first dimensions are plotted!")
    x = val.points[val.vertices,:]
    seriestype --> :shape
    legend --> false
    x[:,1], x[:,2]
end

end
