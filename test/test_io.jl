@testset "File IO" begin
    tempnamefits() do fname
        # create the file first
        FITS(fname, "w") do f
            write(f, ones(2))
        end
        # open as read only
        f = FITS(fname, "r")
        @test summary(f) == "FITS with 1 HDU"
        # writing should throw an error
        @test_throws Exception write(f, ones(2))
        d = f[1]
        @test_throws Exception write(d, ones(2))
        close(f)
        @test summary(f) == "Closed FITS file"

        fname2 = fname*"[12].fits"
        FITS(fname2, "w", extendedparser = false) do f
            write(f, [1,2])
            @test read(f[1]) == [1,2]
        end
        FITSIO.fitswrite(fname2, [1:4;], extendedparser = false)
        @test FITSIO.fitsread(fname2, extendedparser = false) == [1:4;]
    end
end
