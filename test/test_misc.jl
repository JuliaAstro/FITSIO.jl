@test_deprecated FITSIO.libcfitsio_version()

# test we can still access Libcfitsio

@test begin
    using FITSIO.Libcfitsio
    Libcfitsio.libcfitsio_version() isa VersionNumber
end

@testset "0-dim arrays" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            @test_throws ArgumentError write(f, fill(0.0, ()))
            fits_create_empty_img(f.fitsfile)
            CFITSIO.fits_flush_file(f.fitsfile)
            @test isnothing(read(f[1]))
        end
    end
end

@testset "checksum" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            write(f, ones(2,2))
            @test_throws ErrorException FITSIO.verify_checksum(f[1])
            FITSIO.write_checksum(f[1])
            @test FITSIO.verify_checksum(f[1])
            write(f, ones(2,2), checksum=true)
            @test FITSIO.verify_checksum(f[2])
            write(f, Dict("col1" => [1, 2, 3], "col2" => [4, 5, 6]), checksum=true)
            @test FITSIO.verify_checksum(f[3])
            write(f, ["col1", "col2"], [[1,2,3], [4, 5, 6]], checksum=true)
            @test FITSIO.verify_checksum(f[4])
        end
    end
end

@testset "flush" begin
    tempnamefits() do fname
        FITS(fname, "w") do f
            write(f, ones(2,2))
            flush(f)
            f2 = FITS(fname, "r")
            @test length(f2) == 1
            @test read(f2[1]) == ones(2,2)
        end
    end
end
