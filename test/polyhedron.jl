const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include(joinpath(polyhedra_test, "alltests.jl"))

tests = Tuple{String, Function}[]
for n in 2:4
    push!(tests, ("Cross Polytope in $n dimensions", lib->crosspolytopetest(lib, n)))
end

@testset "Polyhedra tests" begin
    runtests(QHull.QHullLib(), tests)
end
