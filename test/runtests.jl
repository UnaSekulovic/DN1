using Test, DN1

@testset "getindex" begin 
    A = RazprsenaMatrika([[3 0 -2 1 0]; [0 3 1 2 0]; [0 0 0 0 0]; [-1 0 0 2 0]; [0 4 0 0 7]])
    @test getindex(A, 1, 3) == -2
    @test_throws BoundsError getindex(A, 0, 0)
    @test getindex(A, 2, 4) == 2
    @test A[3, 4] == 0
    @test A[5, 5] == 7
end

@testset "setindex" begin 
    A = RazprsenaMatrika([[3 0 -2 1 0]; [0 3 1 2 0]; [0 0 0 0 0]; [-1 0 0 2 0]; [0 4 0 0 7]])
    A[1, 3] = 3
    @test A[1, 3] == 3
    A[2, 5] = 3
    @test A[2, 5] == 3
    A[2, 1] = 3
    @test A[2, 1] == 3
    A[4, 2] = 3
    @test A[4, 2] == 3
    A[5, 3] = 3
    @test A[5, 3] == 3
    A[3, 3] = 3
    @test A[3, 3] == 3
end

@testset "firstindex" begin
    A = RazprsenaMatrika([[3 0 -2 1 0]; [0 3 1 2 0]; [0 0 0 0 0]; [-1 0 0 2 0]; [0 4 0 0 7]])
    @test firstindex(A) == (1, 1)
end

@testset "lastindex" begin
    A = RazprsenaMatrika([[3 0 -2 1 0]; [0 3 1 2 0]; [0 0 0 0 0]; [-1 0 0 2 0]; [0 4 0 0 7]])
    @test lastindex(A) == (5, 5)
end

@testset "mnozenje" begin
    A = RazprsenaMatrika([[3 0 -2 1 0]; [0 3 1 2 0]; [0 0 0 0 0]; [-1 0 0 2 0]; [0 4 0 0 7]])
    v = [1, 2, 3, 4, 5]
    @test *(A, v) == [1.0, 17.0, 0.0, 7.0, 43.0]
    v = [0, 0, 0, 0, 0]
    @test *(A, v) == [0.0, 0.0, 0.0, 0.0, 0.0]
    v = [1, -1, 1, -1, 1]
    @test *(A, v) == [0.0, -4.0, 0.0, -3.0, 3.0]
end


@testset "sor" begin
    A = RazprsenaMatrika([[6  -1  0  0  0]; [-1  7  -1 0  0]; [0  -1  8  -1 0]; [0  0  -1  9  -1]; [0  0  0  -1  5]])
    AOsnovna = [[6  -1  0  0  0]; [-1  7  -1 0  0]; [0  -1  8  -1 0]; [0  0  -1  9  -1]; [0  0  0  -1  5]]
    
    x = [1, 1, 1, 1, 1]
    b = AOsnovna * x
    @test AOsnovna \ b ≈ sor(A, b, zeros(size(A.V, 1)), 1)[1]

    x = [1, 2, 3, 4, 5]
    b = AOsnovna * x
    @test AOsnovna \ b ≈ sor(A, b, zeros(size(A.V, 1)), 1.3)[1]

    x = [0, 0, 0, 0, 0]
    b = AOsnovna * x
    @test AOsnovna \ b ≈ sor(A, b, zeros(size(A.V, 1)), 1.4)[1]

end