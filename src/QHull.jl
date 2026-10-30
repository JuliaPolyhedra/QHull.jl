## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module QHull

export Chull, chull

import DirectQhull

mutable struct Chull{T<:Real}
    points::Matrix{T}
    vertices::Vector{Int32}
    simplices::Matrix{Int32}
    facets::Matrix{T}
    area::Float64
    volume::Float64
end

function chull(x::Matrix{T}) where {T<:Real}
    # DirectQhull stores points and facets in columns and uses 1-based indices.
    hull = DirectQhull.ConvexHull(permutedims(Matrix{Float64}(x)))
    vertices = Vector{Int32}(hull.vertices)
    if size(x, 2) > 2
        # Match SciPy's input ordering outside the counterclockwise 2D case.
        sort!(vertices)
    end
    return Chull(
        permutedims(hull.points),
        vertices,
        Matrix{Int32}(permutedims(hull.simplices)),
        permutedims(hull.equations),
        hull.area,
        hull.volume,
    )
end

include("polyhedron.jl")

function Base.show(io::IO, ::MIME"text/plain", ch::Chull)
    println(
        io,
        string(
            "Convex Hull of ",
            size(ch.points, 1),
            " points in ",
            size(ch.points, 2),
            " dimensions",
        ),
    )
    println(io, "Hull segment vertex indices:")
    println(io, ch.vertices)
    println(io, "Points on convex hull in original order:\n")
    println(io, ch.points[sort(ch.vertices[:, 1]), :])
end

end
