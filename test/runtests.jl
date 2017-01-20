using FITSIO
using Base.Test
using Compat

# -----------------------------------------------------------------------------
# Images

# Create a FITS instance and loop over supported types.
fname = tempname() * ".fits"
f = FITS(fname, "w")
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

# test iteration
for hdu in f
    @test size(hdu) == (5, 20)
end
close(f)
if isfile(fname)
    rm(fname)
end

# copy_section()
fname1 = tempname() * ".fits"
f1 = FITS(fname1, "w")
indata = reshape(Float32[1:400;], 20, 20)
write(f1, indata)

fname2 = tempname() * ".fits"
f2 = FITS(fname2, "w")
copy_section(f1[1], f2, 1:10, 1:10)
copy_section(f1[1], f2, 1:10, 1:2:20)
outdata = read(f2[1])
@test outdata == indata[1:10, 1:10]
outdata = read(f2[2])
@test outdata == indata[1:10, 1:2:20]
close(f1)
close(f2)
if isfile(fname1)
    rm(fname1)
end
if isfile(fname2)
    rm(fname2)
end

# -----------------------------------------------------------------------------
# Tables

fname = tempname() * ".fits"
f = FITS(fname, "w")

## Binary table
indata = Dict{Compat.ASCIIString, Array}()
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
indata["vcol1"] = [randstring(j) for j=1:20]  # variable length column
indata["vcol2"] = [collect(1.:j) for j=1.:20.] # variable length

# test writing
write(f, indata; varcols=["vcol1", "vcol2"])

# test reading
for (colname, incol) in indata
    outcol = read(f[2], colname)  # table is in extension 2 (1 = primary hdr)
    @test outcol == incol
    @test eltype(outcol) == eltype(incol)
end


## ASCII tables

indata = Dict{Compat.ASCIIString, Array}()
for (i, T) in enumerate([Int16, Int32, Float32, Float64])
    indata["col$i"] = T[1:20;]
end
i = length(indata) + 1
indata["col$i"] = [randstring(10) for j=1:20]

write(f, indata; hdutype=ASCIITableHDU)

# For ASCII tables, the types don't round trip so we need to define the
# expected output type for each input type.
expected_type = @compat(Dict(Int16=>Int32, Int32=>Int32,
                     Float32=>Float64, Float64=>Float64,
                     Compat.ASCIIString=>Compat.ASCIIString))
for (colname, incol) in indata
    outcol = read(f[3], colname)  # table is in extension 3
    @test outcol == incol
    @test eltype(outcol) == expected_type[eltype(incol)]
end
close(f)
isfile(fname) && rm(fname)

# -----------------------------------------------------------------------------
# FITSHeader

fname = tempname() * ".fits"
f = FITS(fname, "w")

# test that show() works on an empty file and that the beginning and end
# arre what we expect.
io = IOBuffer()
show(io, f)
s = String(take!(io))
@test s[1:6] == "File: "
@test s[end-7:end] == "No HDUs."

inhdr = FITSHeader(["FLTKEY", "INTKEY", "BOOLKEY", "STRKEY", "COMMENT",
                    "HISTORY"],
                   [1.0, 1, true, "string value", nothing, nothing],
                   ["floating point keyword",
                    "",
                    "boolean keyword",
                    "string value",
                    "this is a comment",
                    "this is a history"])

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
s = read_header(f[1], Compat.ASCIIString)
@test s[1:9] == "SIMPLE  ="  # all headers should start with this.
@test length(s) == (9 + length(inhdr)) * 80  # 9 lines = 8 default + "END"

# Test to check that read_header gets the right block even after reading another.
s_reread = read_header(f[1])
s_reread = read_header(f[2])
s_reread = read_header(f[1], Compat.ASCIIString)
@assert s == s_reread

# Read single keywords
@test read_key(f[1], 9) == ("FLTKEY", 2.0, "floating point keyword")
@test read_key(f[1], "FLTKEY") == (2.0, "floating point keyword")

# Test that show() works and that the beginning of output is what we expect.
io = IOBuffer()
show(io, f)
s = String(take!(io))
@test s[1:6] == "File: "

# Clean up from last test.
close(f)
if isfile(fname)
    rm(fname)
end

# -----------------------------------------------------------------------------
# parsing non-standard keyword records

# `create_test_file` : Create a simple FITS file for testing, with the
# given header string added after the required keywords. The length of
# `header` must be a multiple of 80.  The purpose of creating such
# files is to test the parsing of non-standard FITS keyword records
# (non-standard files can't be created with cfitsio).

function create_test_file(fname::AbstractString, header::Compat.ASCIIString)
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

# Create a test file with a few non-standard keyword records
fname = tempname() * ".fits"
header = "        Warning: CROTA2 is inaccurate due to considerable skew                  SKEW    =  1.9511305786508E+00,  1.9234924037208E+00 /Measure of skew           "
create_test_file(fname, header)

# check that we can read the header (and data).
f = FITS(fname)
hdr = read_header(f[1])
data = read(f[1])
close(f)

# test that we can write it back out.
fname2 = tempname() * ".fits"
f2 = FITS(fname2, "w")
write(f2, data; header=hdr)
close(f2)

# clean up.
if isfile(fname)
    rm(fname)
end
if isfile(fname2)
    rm(fname2)
end

# -----------------------------------------------------------------------------
# Miscellaneous

# test that this function works and returns the right type.
@test typeof(FITSIO.libcfitsio_version()) === VersionNumber

println("All tests passed.")
