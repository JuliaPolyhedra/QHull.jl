using Test
import QHull
import Polyhedra
using LinearAlgebra

@testset "chull" begin
    pts = [
        -1.0 0
        1 1
        3 0
        3 3
        2 2
        0 3
        -1 2
    ]

    hull = QHull.chull(pts)
    @test hull.vertices == [1, 3, 4, 6, 7]
    @test size(hull.points) == size(pts)
    @test hull.simplices == Int32[3 1; 4 3; 6 4; 7 1; 7 6]

    # multi-dim
    x = randn(1000, 5)
    hull = QHull.chull(x)
end

# Check the row-based public API against known geometry, including integer and
# Float32 inputs which the native backend itself does not accept.
@testset "Native convex hull" begin
    for T in (Int, Float32, Float64)
        points = T[0 0; 2 0; 2 2; 0 2; 1 1; 0 0]
        hull = QHull.chull(points)
        @test hull isa QHull.Chull{Float64}
        @test hull.points == points
        @test hull.points !== points
        @test hull.vertices == Int32[1, 2, 3, 4]
        @test eltype(hull.simplices) == Int32
        @test size(hull.simplices) == (4, 2)
        @test hull.area ≈ 8
        @test hull.volume ≈ 4
        @test all(points * hull.facets[:, 1:2]' .+ hull.facets[:, 3]' .<= 1e-12)
        for i in axes(hull.simplices, 1)
            @test points[hull.simplices[i, :], :] * hull.facets[i, 1:2] .+
                  hull.facets[i, 3] ≈ zeros(2) atol=1e-12
        end
    end

    @testset "Simplex in $d dimensions" for d in (3, 5)
        points = vcat(zeros(1, d), Matrix{Float64}(I, d, d), fill(1 / (2d), 1, d))
        hull = QHull.chull(points)
        @test hull.points == points
        @test hull.vertices == Int32.(1:(d+1))
        @test size(hull.simplices) == (d+1, d)
        @test size(hull.facets) == (d+1, d+1)
        @test hull.volume ≈ 1 / factorial(d)
        @test hull.area ≈ (d + sqrt(d)) / factorial(d-1)
        @test all(points * hull.facets[:, 1:d]' .+ hull.facets[:, d+1]' .<= 1e-12)
        for i in axes(hull.simplices, 1)
            @test points[hull.simplices[i, :], :] * hull.facets[i, 1:d] .+
                  hull.facets[i, d+1] ≈ zeros(d) atol=1e-12
        end
    end
end

@testset "Native Polyhedra conversions" begin
    # Use explicit rows for the eight cube vertices.
    points = [collect(p) for p in Iterators.product((-1.0, 1.0), (-1.0, 1.0), (-1.0, 1.0))]
    points = reduce(vcat, permutedims.(vec(points)))
    p = Polyhedra.polyhedron(Polyhedra.vrep(points), QHull.Library())
    @test Polyhedra.volume(p) ≈ 8
    @test Polyhedra.nhalfspaces(Polyhedra.hrep(p)) == 6
    @test Polyhedra.npoints(Polyhedra.vrep(p)) == 8
    h = Polyhedra.hrep(vcat(Matrix{Float64}(I, 3, 3), -Matrix{Float64}(I, 3, 3)), ones(6))
    p = Polyhedra.polyhedron(h, QHull.Library())
    @test Polyhedra.npoints(Polyhedra.vrep(p)) == 8
    @test Polyhedra.volume(p) ≈ 8
end
