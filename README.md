QHull
=====
[![Build Status](https://travis-ci.org/davidavdav/QHull.jl.svg)](https://travis-ci.org/davidavdav/QHull.jl)


A Julia wrapper around a PyCall wrapper around `scipy.spatial.ConvexHull`, which uses the qhull Convex Hull library.

The qhull library for computing the convex hull of data points seems to be the standard and very widely used.

This module is a quick wrapper around a Python wrapper around the library, as suggested by [Miles Lubin](https://groups.google.com/d/topic/julia-users/e9m8t5W3TVs/discussion).

Synopsis
--------

```julia

using QHull

p = rand(10,2)
ch = chull(p)
ch.points         # original points
ch.vertices       # indices to line segments forming the convex hull
ch.simplices      # the simplexes forming the convex hull
show(ch)
```
