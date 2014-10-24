## CHull.jl 
## (c) 2013 David Al van Leeuwen
## A Julia wrapper around a python wrapper around the qhull Convex Hull library

## This code is licensed under the GNU General Public License, version 2
## See the file LICENSE in this distribution

module CHull

export Chull, chull, display, show

using PyCall
@pyimport pyhull
@pyimport pyhull.convex_hull as convex_hull

type Chull{T<:Real}
    points::Array{T}
    vertices::Array{Int}
end

function chull{T<:Real}(x::Array{T})
    py = convex_hull.ConvexHull(x)
    points = convert(Array,py["points"])[1]
    vertices = convert(Array{Int},py["vertices"]) + 1
    Chull(points, vertices)
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

end
