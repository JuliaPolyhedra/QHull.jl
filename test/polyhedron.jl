using Polyhedra
const polyhedra_test = joinpath(dirname(dirname(pathof(Polyhedra))), "test")

using JuMP
import GLPK
# Need `"presolve" => GLPK.ON` for `detect_new_linearities`, see
# https://travis-ci.org/github/JuliaPolyhedra/Polyhedra.jl/jobs/691916637#L486
lp_solver = optimizer_with_attributes(GLPK.Optimizer, "presolve" => GLPK.GLP_ON, MOI.Silent() => true)
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
    "orthantdecompose", "largedecompose", "recipe", "issue224",
    # 1D not supported:
    "vhypercubetest1c", "vhypercubetest1u"
]

polyhedratest(LIB, exclude)
