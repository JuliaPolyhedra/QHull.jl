const polyhedra_test = joinpath(Pkg.dir("Polyhedra"), "test")

include("solvers.jl")

if !isempty(lp_solvers)
    lpsolver = first(lp_solvers)
end

include(joinpath(polyhedra_test, "alltests.jl"))

if !isempty(lp_solvers)
    tests = Tuple{String, Function}[]
    for n in 2:2
        push!(tests, ("Hypercube in $n dimensions", lib->hypercubetest(lib, n)))
        push!(tests, ("Simplex with the origin in $n dimensions", lib->simplexorigtest(lib, n)))
        push!(tests, ("Cross Polytope in $n dimensions", lib->crosspolytopetest(lib, n)))
    end
    push!(tests, ("Doc example", doctest))

    @testset "Polyhedra tests" begin
        runtests(QHull.QHullLibrary(lpsolver), tests)
    end
end
