using Test
import QHull

include("chull.jl")

@testset "Polyhedra tests" begin
    include("polyhedron.jl")
end
