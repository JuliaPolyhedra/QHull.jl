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
    area::T
    equations::Matrix{T}
    ndim::Int
    max_bound::Vector{T}
    min_bound::Vector{T}
    neighbours::Matrix{T}
    npoints::Int
    nsimplex::Int
    points::Matrix{T}
    simplices::Vector{Vector{Int}}
    vertices::Vector{Int}
    volume::T
end


function chull{T<:Real}(x::Matrix{T})

    hull = spatial[:ConvexHull](x)

    ch = Chull([hull[field] for field in fieldnames(Chull)]...)

    # fix base-0 / base-1 difference
    ch.vertices += 1
    ch.simplices += 1
    return ch
end

include("polyhedron.jl")

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
