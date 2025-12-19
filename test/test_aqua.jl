using Aqua, FITSIO

@testset "project quality" begin
    Aqua.test_all(FITSIO)
end
