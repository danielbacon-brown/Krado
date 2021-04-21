@testset "NestedArrays" begin

    A = [1.0 4.0; 5.0 9.0]
    B = [7.0 3.0;-4.0 0.0]
    C = [2.0 -5.0;1.0 8.0]
    D = [5.0 3.0; 4.0 7.0]

    # Test creation of nested matrix:
    emptyArr = zeros(Float64, 0,0)
    ABCD_split_man = fill(emptyArr,2,2)
    ABCD_split_man[1,1] = A
    ABCD_split_man[1,2] = B
    ABCD_split_man[2,1] = C
    ABCD_split_man[2,2] = D
    ABCD_split = [ [A] [B]; [C] [D] ]
    @test all( isapprox( ABCD_split, ABCD_split_man, atol=1e-3) )

    # Test conversion of split array into combined array
    ABCD_recombined = combineArray(ABCD_split)
    ABCD_recombined_man = [ ABCD_split[1,1] ABCD_split[1,2]; ABCD_split[2,1] ABCD_split[2,2] ]
    @test all( isapprox( ABCD_recombined, ABCD_recombined_man, atol=1e-3 ) )

    # Test conversion of combined array into split array
    ABCD_resplit = splitArray(ABCD_recombined)
    ABCD_split = [ [A] [B]; [C] [D] ]
    @test all( isapprox( ABCD_resplit, ABCD_split, atol=1e-3 ) )

    # Test multiplication of 2x2 nested arrays
    ABCD = [A B; C D]
    DCBA = [D C; B A]
    ABCD_nest = splitArray(ABCD)
    DCBA_nest = splitArray(DCBA)
    @test all( isapprox( ABCD*DCBA, combineArray(ABCD_nest*DCBA_nest) ))

    # Test RedhefferStarProduct:
    RedhefferStarProd = RedhefferStarProduct(ABCD_nest, DCBA_nest)
    refProd = [134.043 115.13 -0.944444 5.86111; -45.087 2.73913 0.777778 -5.44444; -1.08333 8.83333 134.043 115.13; 0.916667 -8.41667 -45.087 2.73913]
    @test all( isapprox( combineArray(RedhefferStarProd), refProd, atol=1e-3))

end;
