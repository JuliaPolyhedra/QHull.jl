using Polyhedra
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(polyhedra_test, "solvers.jl"))
include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))

exclude = [
    "jumpsimplex", "ex1", "infeasible", "nonfulldimensional",
    "simplex", "permutahedron", "board", "issue48", "empty",
    "orthantdecompose", "largedecompose", "recipe",
    # 1D not supported:
    "vhypercubetest1c", "vhypercubetest1u"
]

polyhedratest(QHull.Library(lp_solver), exclude)
