# -----------------------------------------------------------------------------
# parsing non-standard keyword records

# `create_test_file` : Create a simple FITS file for testing, with the
# given header string added after the required keywords. The length of
# `header` must be a multiple of 80.  The purpose of creating such
# files is to test the parsing of non-standard FITS keyword records
# (non-standard files can't be created with cfitsio).

function create_test_file(fname::AbstractString, header::String)
    if length(header) % 80 != 0
        error("length of header must be multiple of 80")
    end

    open(fname, "w") do f

        stdhdr = "SIMPLE  =                    T / file does conform to FITS standard             BITPIX  =                  -64 / number of bits per data pixel                  NAXIS   =                    2 / number of data axes                            NAXIS1  =                   10 / length of data axis 1                          NAXIS2  =                   10 / length of data axis 2                          EXTEND  =                    T / FITS dataset may contain extensions            "
        endline = "END                                                                             "
        data = fill(0., (10, 10))  # 10x10 array of big-endian Float64 zeros

        # write header
        write(f, stdhdr)
        write(f, header)
        write(f, endline)

        # add padding
        block_position = (length(stdhdr) + length(header) + length(endline)) % 2880
        padding = (block_position == 0) ? 0 : 2880 - block_position
        write(f, " "^padding)

        # write data
        write(f, data)

        # add padding
        block_position = sizeof(data) % 2880
        padding = (block_position == 0) ? 0 : 2880 - block_position
        write(f, fill(0x00, (padding,)))
    end
end

@testset "Non-standard keywords" begin
    mktempdir() do dir
        fname = randomfilename(dir)
        # Create a test file with a few non-standard keyword records
        header = "        Warning: CROTA2 is inaccurate due to considerable skew                  SKEW    =  1.9511305786508E+00,  1.9234924037208E+00 /Measure of skew           "
        create_test_file(fname, header)

        # check that we can read the header (and data).
        FITS(fname) do f
            hdr = read_header(f[1])
            @test hdr["SKEW"] == "1.9511305786508E+00,"
            @test get_comment(hdr, "SKEW") == "1.9234924037208E+00 /Measure of skew"
            data = read(f[1])
            @test data == zeros(10, 10)
        end

        # Test that we can read is at a memory backed file
        FITS(read(fname)) do f
            hdr = read_header(f[1])
            data = read(f[1])
            @test data == zeros(10, 10)
            # test that we can write it back out.
            FITS(randomfilename(dir), "w") do f2
                write(f2, data; header=hdr)
            end
        end
    end
end
