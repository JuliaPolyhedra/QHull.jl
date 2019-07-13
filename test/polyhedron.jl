using Polyhedra
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

include(joinpath(polyhedra_test, "solvers.jl"))

const LIB = QHull.Library(lp_solver)

@testset "Volume" begin
    h = HalfSpace([1, 0, 0], 1) ∩ HalfSpace([0, 1, 0], 1) ∩ HalfSpace([0, 0, 1], 1) ∩
        HalfSpace([-1, 0, 1], 1) ∩ HalfSpace([0, -1, 0], 1) ∩ HalfSpace([0, 0, -1], 1)
    p = polyhedron(h, LIB)
    @test volume(p) ≈ 8
end

include(joinpath(polyhedra_test, "utils.jl"))
include(joinpath(polyhedra_test, "polyhedra.jl"))

exclude = [
    "jumpsimplex", "ex1", "infeasible", "nonfulldimensional",
    "simplex", "permutahedron", "board", "issue48", "empty",
    "orthantdecompose", "largedecompose", "recipe",
    # 1D not supported:
    "vhypercubetest1c", "vhypercubetest1u"
]

polyhedratest(LIB, exclude)
