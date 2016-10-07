#!/usr/bin/env julia
include("src/CHull.jl")
using CHull

cd("test")
include("test/runtests.jl")
