FITSIO.jl
=========

Flexible Image Transport System (FITS) support for Julia


[![Build Status](https://github.com/juliaastro/FITSIO.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/FITSIO.jl/actions)
[![PkgEval](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/F/FITSIO.svg)](https://juliaci.github.io/NanosoldierReports/pkgeval_badges/report.html)
[![codecov](https://codecov.io/gh/juliaastro/fitsio.jl/branch/master/graph/badge.svg?token=SA9EG0z8pt)](https://codecov.io/gh/juliaastro/fitsio.jl)


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://juliaastro.github.io/FITSIO.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](http://juliaastro.github.io/FITSIO.jl/dev/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Note:** The `Libcfitsio` submodule has been moved to [CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl) and will be deprecated in a future release.


## Usage

For more in-depth usage and examples, see [the documentation](http://juliaastro.github.io/FITSIO.jl/stable/).
Here, we provide an example where we read a sample fits file [provided by NASA](https://fits.gsfc.nasa.gov/fits_samples.html).

```julia
julia> using FITSIO, Downloads

julia> fname = Downloads.download("https://fits.gsfc.nasa.gov/samples/FOSy19g0309t_c2f.fits")

julia> f = FITS(fname, "r")
File: docs/src/FOSy19g0309t_c2f.fits
Mode: "r" (read-only)
HDUs: Num  Name               Type
      1                       Image
      2    y19g0309t.c2h.tab  ASCIITable

julia> f[1]
File: docs/src/FOSy19g0309t_c2f.fits
HDU: 1
Mode: read-only
Type: Image
Datatype: Float32
Datasize: (2064, 2)

# read an image from disk
julia> data = read(f[1]);

# read just a subset of image
julia> read(f[1], 1:4, :)
4×2 Matrix{Float32}:
 2.48511f-15  1.36115f-15
 1.56953f-15  1.10982f-15
 0.0          0.0
 1.12148f-15  9.71231f-16

# Get info about a table
julia> f[2]
File: docs/src/FOSy19g0309t_c2f.fits
HDU: 2 (name=y19g0309t.c2h.tab)
Type: ASCIITable
Rows: 2
Columns: Name      Type     TFORM
         CRVAL1    Float64  D25.16
         CRPIX1    Float64  E15.7
         CD1_1     Float64  E15.7
         DATAMIN   Float64  E15.7
         DATAMAX   Float64  E15.7
         RA_APER   Float64  D25.16
         DEC_APER  Float64  D25.16
         FILLCNT   Int32    I11
         ERRCNT    Int32    I11
         FPKTTIME  Float64  D25.16
         LPKTTIME  Float64  D25.16
         CTYPE1    String   A8
         APER_POS  String   A8
         PASS_DIR  Int32    I11
         YPOS      Float64  E15.7
         YTYPE     String   A4
         EXPOSURE  Float64  E15.7
         X_OFFSET  Float64  E15.7
         Y_OFFSET  Float64  E15.7

# Read a column from the table:
julia> data = read(f[2], "CD1_1")
2-element Vector{Float64}:
 1.0
 1.0

# Read the entire header into memory
julia> header = read_header(f[1]);

julia> header["NAXIS1"]  # get value by keyword
2064

# Read single keys into memory
julia> read_key(f[1], 4) # by position
("NAXIS1", 2064, "")

julia> read_key(f[1], "NAXIS1") # by keyword
(2064, "")

# write data to file
julia> FITS("new_file.fits", "w") do f
           write(f, data; header)
       end
```
