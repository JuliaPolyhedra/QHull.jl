QHull
=====
[![Build Status](https://travis-ci.org/JuliaPolyhedra/QHull.jl.svg)](https://travis-ci.org/JuliaPolyhedra/QHull.jl)


A Julia wrapper around a PyCall wrapper around `scipy.spatial.ConvexHull`, which uses the qhull Convex Hull library.
It implements the Polyhedral Computation library interface of [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl).

The qhull library for computing the convex hull of data points seems to be the standard and very widely used.

This module is a quick wrapper around a Python wrapper around the library, as suggested by [Miles Lubin](https://groups.google.com/d/topic/julia-users/e9m8t5W3TVs/discussion).

Synopsis
--------

Low-level interface:
```julia
using QHull

p = rand(10,2)
ch = chull(p)
ch.points         # original points
ch.vertices       # indices to line segments forming the convex hull
ch.simplices      # the simplexes forming the convex hull
show(ch)
```

Using [Polyhedra.jl](https://github.com/JuliaPolyhedra/Polyhedra.jl):
```julia
using Polyhedra
# Constructs a V-representation of 10 random points in 2 dimension
v = vrep(rand(10, 2))

using QHull
# Constructs a polyhedon from this V-representation with the QHull library
p = polyhedron(v, QHull.Library())
# Removing redundant points, i.e. points which are in the interior of the convex hull
removevredundancy!(p)
# Show remaining points, i.e. the non-redundant ones
removevredundancy!(p)
@show vrep(p)
# Show the H-representation, the facets describing the polytope
@show hrep(p)
```
