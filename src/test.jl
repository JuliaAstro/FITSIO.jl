using FITSIO
using Base.Test

# images of various types
fname = tempname() * ".fits"
f = FITS(fname, "w")
for T in [Uint8, Int8, Uint16, Int16, Uint32, Int32, Int64,
          Float32, Float64]
    indata = reshape(T[1:100], 5, 20)
    write(f, indata)
    outdata = read(f[end])
    @test indata == outdata
    @test eltype(indata) == eltype(outdata)
    end
end
close(f)
if isfile(fname)
    rm(fname)
end
