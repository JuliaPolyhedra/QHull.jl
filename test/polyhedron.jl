using Polyhedra
const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "solvers.jl"))
include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))

@testset "Polyhedra tests" begin
    @show lp_solver
    polyhedratest(QHull.QHullLibrary(lp_solver), ["jumpsimplex", "ex1", "infeasible", "nonfulldimensional",
                                                  "doc", "simplex", "permutahedron", "board", "issue48", "empty",
                                                  "orthantdecompose", "largedecompose"])
end
