using Polyhedra
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(polyhedra_test, "solvers.jl"))
include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))

@testset "Polyhedra tests" begin
    polyhedratest(QHull.Library(lp_solver), ["jumpsimplex", "ex1", "infeasible", "nonfulldimensional",
                                             "simplex", "permutahedron", "board", "issue48", "empty",
                                             "orthantdecompose", "largedecompose", "recipe"])
end
