## CHull.jl
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module CHull

export Chull, chull, display, show

using PyCall
@pyimport scipy.spatial as spatial

type Chull{T<:Real}
    points::Array{Array{T,1},1}
    vertices::Array{Int,1}
    simplices::Array{Any,1}
end

function chull{T<:Real}(x::Array{T})
    py = spatial.ConvexHull(x)
    points = convert(Array{Array{T,1},1},py["points"])
    vertices = convert(Array{Int},py["vertices"]) + 1
    simplices = convert(Array{Array{Int,1},1},py["simplices"]) + 1
    Chull(points, vertices, simplices)
end

import Base.show
## I don't seem to be able to print newlines in the display()
#function display(t::TextDisplay, ch::Chull)
#    display(t, string("Convex Hull of ", size(ch.points,1), " points in ", size(ch.points,2), "dimensions"))
#    display(t, ch.vertices)
#    display(t, ch.points[sort(ch.vertices[:,1]),:])
#end

## should I use print statements in show()
function show(io::IO, ch::Chull)
    println(io, string("Convex Hull of ", size(ch.points,1), " points in ", size(ch.points,2), " dimensions"))
    println(io, "Hull segment vertex indices:")
    println(io, ch.vertices)
    println(io, "Points on convex hull in original order:\n")
    println(io, ch.points[sort(ch.vertices[:,1]),:])
end

using RecipesBase
@recipe function f{T<:Chull}(val::T)
    length(val.points[1]) && warning("Only the two first dimensions are plotted!")
    x = [x[1] for x in val.points]
    y = [x[2] for x in val.points]
    seriestype --> :shape
    legend --> false
    x[val.vertices], y[val.vertices]
end

end
