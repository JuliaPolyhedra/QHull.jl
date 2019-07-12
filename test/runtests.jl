using Test
import QHull

@testset "chull" begin
    pts = [-1.0 0;
           1 1;
           3 0;
           3 3;
           2 2;
           0 3;
           -1 2]

    hull = QHull.chull(pts)
    @test hull.vertices == [1, 3, 4, 6, 7]
    @test size(hull.points) == size(pts)
    @test hull.simplices == Array{Int,1}[[3,1], [4,3], [6,4], [7,1], [7,6]]

    # multi-dim
    x = randn(1000, 5)
    hull = QHull.chull(x)
end

@testset "Polyhedra tests" begin
    include("polyhedron.jl")
end
