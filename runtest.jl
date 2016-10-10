#!/usr/bin/env julia
include("src/QHull.jl")
using QHull

cd("test")
include("test/runtests.jl")
