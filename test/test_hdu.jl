@testset "Write data to an existing image HDU" begin
    @testset "Overwrite an image" begin
        # Create a file with two images
        tempnamefits() do fname

            FITS(fname, "w") do f
                write(f, [[1 2 3]; [4 5 6]])
                write(f, [[7 8 9]; [17 52 10]])
            end

            # Open the file in read-write mode
            FITS(fname, "r+") do f
                @test read(f[1]) == [[1 2 3]; [4 5 6]]
                @test read(f[2]) == [[7 8 9]; [17 52 10]]

                # Write data into the first image HDU
                write(f[1], [[11 12 13]; [14 15 16]])

                @test read(f[1]) == [[11 12 13]; [14 15 16]]  # First HDU is updated
                @test read(f[2]) == [[7 8 9]; [17 52 10]]  # Second HDU is untouched
            end
        end
    end

    @testset "Show size mismatch error" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, [[1 2 3 4]; [5 6 7 8]])
                write(f, [[1 2 3]; [4 5 6]])
            end

            FITS(fname, "r+") do f
                second_hdu = f[2]  # Get the second HDU

                # Make the first HDU to be the current HDU
                # Needed to check if `write` method compares data size with the given HDU
                # and not with the current HDU
                fits_movabs_hdu(f.fitsfile, 1)

                # Attempt to write data array of incorrect size to the second HDU
                exception = ErrorException("size of HDU (2, 3) not equal to size of data (2, 4).")
                @test_throws exception write(second_hdu, [[11 12 13 14]; [15 16 17 18]])
            end
        end
    end

    @testset "write with different eltype" begin
        tempnamefits() do fname
            FITS(fname, "w") do f
                write(f, ones(2,2)) # create an image with eltype Float64
                write(f[1], ones(Int, 2, 2) .* 2) # write Int data to the same HDU
                @test read(f[1]) == ones(2,2) .* 2 # check that data is written correctly
            end
        end
    end

    @testset "fitswrite" begin
        tempnamefits() do fname
            a = ones(3,3)
            FITSIO.fitswrite(fname, a)
            FITS(fname, "r") do f
                @test read(f[1]) == a
            end
        end
    end
end
