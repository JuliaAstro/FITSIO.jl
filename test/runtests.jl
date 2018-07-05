using FITSIO
using Test

@testset "Images" begin
    # Create a FITS instance and loop over supported types.
    fname = tempname() * ".fits"
    FITS(fname, "w") do f
        for T in [UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64,
                  Float32, Float64]
            indata = reshape(T[1:100;], 5, 20)

            # Test writing the data to a new extension
            write(f, indata)

            # test reading the full array
            outdata = read(f[end])
            @test indata == outdata
            @test eltype(indata) == eltype(outdata)

            # test reading subsets of the array
            @test read(f[end], :, :) == indata
            @test read(f[end], 4, 1:10) == indata[4, 1:10]  # 2-d array
            @test read(f[end], :, 4) == indata[:, 4]  # 1-d array
            @test read(f[end], 2, 3) == indata[2, 3]  # scalar
            @test read(f[end], :, 1:2:10) == indata[:, 1:2:10]
            @test read(f[end], 1:3, :) == indata[1:3, :]

            # test expected errors
            @test_throws DimensionMismatch read(f[end], :)
            @test_throws DimensionMismatch read(f[end], :, :, 1)
            @test_throws BoundsError read(f[end], 1:6, :)
            @test_throws BoundsError read(f[end], 1, 0)

        end

        @test_throws ErrorException f[100]

        # Test representation
        @test repr(f)[end-17:end] == "9          Image  "
        @test repr(f[1])[1:6] == "File: "

        # test iteration
        for hdu in f
            @test size(hdu) == (5, 20)
        end
    end
    rm(fname, force=true)

    # copy_section()
    fname1 = tempname() * ".fits"
    fname2 = tempname() * ".fits"
    try
        f1 = FITS(fname1, "w")
        indata = reshape(Float32[1:400;], 20, 20)
        write(f1, indata)

        f2 = FITS(fname2, "w")
        copy_section(f1[1], f2, 1:10, 1:10)
        copy_section(f1[1], f2, 1:10, 1:2:20)
        outdata = read(f2[1])
        @test outdata == indata[1:10, 1:10]
        outdata = read(f2[2])
        @test outdata == indata[1:10, 1:2:20]
        close(f1)
        close(f2)
    finally
        rm(fname1, force=true)
        rm(fname2, force=true)
    end
end

@testset "Tables" begin
    fname = tempname() * ".fits"
    f = FITS(fname, "w")

    ## Binary table
    indata = Dict{String, Array}()
    for (i, T) in enumerate([UInt8, Int8, UInt16, Int16, UInt32, Int32, Int64,
                             Float32, Float64, Complex64, Complex128])
        indata["col$i"] = T[1:20;]
    end
    i = length(indata) + 1
    indata["col$i"] = [randstring(10) for j=1:20]  # ASCIIString column
    i += 1
    indata["col$i"] = ones(Bool, 20)  # Bool column
    i += 1
    indata["col$i"] = reshape([1:40;], (2, 20))  # vector Int64 column
    i += 1
    indata["col$i"] = [randstring(5) for j=1:2, k=1:20]  # vector ASCIIString col
    indata["vcol"] = [randstring(j) for j=1:20]  # variable length column
    indata["VCOL"] = [collect(1.:j) for j=1.:20.] # variable length

    # test writing
    write(f, indata; varcols=["vcol", "VCOL"])

    # test reading
    colnames = FITSIO.colnames(f[2])
    for (colname, incol) in indata
        outcol = read(f[2], colname)  # table is in extension 2 (1 = primary hdr)
        @test outcol == incol
        @test eltype(outcol) == eltype(incol)
        @test colname in colnames
    end

    @test_throws ErrorException read(f[2], "vcol", case_sensitive=false)

    # Test representation
    @test repr(f[2])[end-38:end] == "\n\n         (*) = variable-length column"

    ## ASCII tables

    indata = Dict{String, Array}()
    for (i, T) in enumerate([Int16, Int32, Float32, Float64])
        indata["col$i"] = T[1:20;]
    end
    i = length(indata) + 1
    indata["col$i"] = [randstring(10) for j=1:20]

    write(f, indata; hdutype=ASCIITableHDU)

    # For ASCII tables, the types don't round trip so we need to define the
    # expected output type for each input type.
    expected_type = Dict(Int16=>Int32, Int32=>Int32,
                         Float32=>Float64, Float64=>Float64,
                         String=>String)
    colnames = FITSIO.colnames(f[3])
    for (colname, incol) in indata
        outcol = read(f[3], colname)  # table is in extension 3
        @test outcol == incol
        @test eltype(outcol) == expected_type[eltype(incol)]
        @test colname in colnames
    end
    # test show/repr on ASCIITableHDU by checking that a couple lines are what we expect
    lines = split(repr(f[3]), "\n")
    @test lines[4] == "Rows: 20"
    @test lines[6] == "         col3  Float64  E26.17  "
    close(f)
    rm(fname, force=true)
end

@testset "FITSHeader" begin
    fname = tempname() * ".fits"
    f = FITS(fname, "w")

    # test that show() works on an empty file and that the beginning and end
    # arre what we expect.
    s = repr(f)
    @test s[1:6] == "File: "
    @test s[end-7:end] == "No HDUs."

    @test_throws ErrorException FITSHeader(["KEY1"], [1, 2, 3], ["comment 1", "comment 2"])

    inhdr = FITSHeader(["FLTKEY", "INTKEY", "BOOLKEY", "STRKEY", "COMMENT",
                        "HISTORY"],
                       [1.0, 1, true, "string value", nothing, nothing],
                       ["floating point keyword",
                        "",
                        "boolean keyword",
                        "string value",
                        "this is a comment",
                        "this is a history"])

    @test repr(inhdr) == """
FLTKEY  =                  1.0 / floating point keyword
INTKEY  =                    1
BOOLKEY =                    T / boolean keyword
STRKEY  =       'string value' / string value
COMMENT this is a comment
HISTORY this is a history"""

    inhdr["INTKEY"] = 2  # test setting by key
    inhdr[1] = 2.0  # test settting by index
    set_comment!(inhdr, "INTKEY", "integer keyword") # test setting a comment

    indata = reshape(Float32[1:100;], 5, 20)
    write(f, indata; header=inhdr)

    # Write a second block.
    inhdr2 = deepcopy(inhdr)
    inhdr2["INTKEY"] = 3 # Set it to a different value.
    write(f, indata; header=inhdr2)

    outhdr = read_header(f[1])
    @test outhdr["FLTKEY"] === 2.0
    @test outhdr["INTKEY"] === 2
    @test outhdr["BOOLKEY"] === true
    @test outhdr["STRKEY"] == "string value"
    @test get_comment(outhdr, 13) == "this is a comment"
    @test get_comment(outhdr, 14) == "this is a history"
    @test length(outhdr) == 14
    @test haskey(outhdr, "FLTKEY")

    # Read entire header as a single string
    s = read_header(f[1], String)
    @test s[1:9] == "SIMPLE  ="  # all headers should start with this.
    @test length(s) == (9 + length(inhdr)) * 80  # 9 lines = 8 default + "END"

    # Test to check that read_header gets the right block even after reading another.
    s_reread = read_header(f[1])
    s_reread = read_header(f[2])
    s_reread = read_header(f[1], String)
    @test s == s_reread

    # update an existing keyword, and read it directly
    write_key(f[1], "FLTKEY", 3.0)
    @test read_key(f[1], 9) == ("FLTKEY", 3.0, "floating point keyword")
    @test read_key(f[1], "FLTKEY") == (3.0, "floating point keyword")

    # Test appending a keyword, then modifying a keyword of different
    # values with write_key()
    for value in [1.0, "string value", 42, false, nothing]
        write_key(f[1], "NEWKEY", value, "new key comment")
        @test read_key(f[1], "NEWKEY") == (value, "new key comment")
        @test read_key(f[1], 15) == ("NEWKEY", value, "new key comment")
    end

    # Test that show() works and that the beginning of output is what we expect.
    @test repr(f)[1:6] == "File: "

    # Clean up from last test.
    close(f)
    rm(fname, force=true)
end

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

    f = open(fname, "w")

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

    close(f)
end

@testset "Non-standard keywords" begin
    fname = tempname() * ".fits"
    fname2 = tempname() * ".fits"
    try
        # Create a test file with a few non-standard keyword records
        header = "        Warning: CROTA2 is inaccurate due to considerable skew                  SKEW    =  1.9511305786508E+00,  1.9234924037208E+00 /Measure of skew           "
        create_test_file(fname, header)

        # check that we can read the header (and data).
        f = FITS(fname)
        hdr = read_header(f[1])
        @test hdr["SKEW"] == "1.9511305786508E+00,"
        @test get_comment(hdr, "SKEW") == "1.9234924037208E+00 /Measure of skew"
        data = read(f[1])
        @test data == zeros(10, 10)
        close(f)

        # Test that we can read is at a memory backed file
        f = FITS(read(fname))
        hdr = read_header(f[1])
        data = read(f[1])
        @test data == zeros(10, 10)
        close(f)

        # test that we can write it back out.
        f2 = FITS(fname2, "w")
        write(f2, data; header=hdr)
        close(f2)
    finally
        rm(fname, force=true)
        rm(fname2, force=true)
    end
end

# -----------------------------------------------------------------------------

@testset "Miscellaneous" begin
    # test that this function works and returns the right type.
    @test typeof(FITSIO.libcfitsio_version()) === VersionNumber
    # test it parses a number as intended.
    @test FITSIO.libcfitsio_version(3.341)  === VersionNumber(3, 34, 1)
    @test FITSIO.libcfitsio_version(3.41f0) === VersionNumber(3, 41, 0)
end

@testset "Libcfitsio" begin
    include("libcfitsio.jl")
end
