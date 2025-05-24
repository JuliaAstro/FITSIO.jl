# FITSIO.jl

[![GitHub](https://img.shields.io/badge/Code-GitHub-black.svg)](https://github.com/JuliaAstro/FITSIO.jl)
[![Build Status](https://github.com/juliaastro/FITSIO.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/FITSIO.jl/actions)
[![Coverage Status](http://img.shields.io/coveralls/JuliaAstro/FITSIO.jl.svg?style=flat-square)](https://coveralls.io/r/JuliaAstro/FITSIO.jl?branch=master)

A [Julia](http://julialang.org) package for reading and writing
Flexible Image Transport System (FITS) files, based on the
[cfitsio](http://heasarc.gsfc.nasa.gov/fitsio/) library.

The interface is inspired by Erin Sheldon's
[fitsio](https://github.com/esheldon/fitsio) Python package.

!!! warning
    The `Libcfitsio` submodule has been moved to [CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl) and will be deprecated in a future release.

## Installation

`FITSIO.jl` can be installed using the built-in package manager

```julia-repl
pkg> add FITSIO
```

## Usage

## Reading a file

We read [a FITS file](https://fits.gsfc.nasa.gov/samples/FOSy19g0309t_c2f.fits) from [NASA's list of sample files](https://fits.gsfc.nasa.gov/fits_samples.html).

```@meta
DocTestFilters = r"File: [a-zA-Z0-9/\.]*docs/src/FOSy19g0309t_c2f.fits"
```

To open an existing file for reading:
```jldoctest nasa_sample
julia> using FITSIO

julia> filepath = joinpath(dirname(pathof(FITSIO)), "../docs/src/FOSy19g0309t_c2f.fits");

julia> f = FITS(filepath, "r")
File: docs/src/FOSy19g0309t_c2f.fits
Mode: "r" (read-only)
HDUs: Num  Name               Type
      1                       Image
      2    y19g0309t.c2h.tab  ASCIITable
```

At the REPL, information about the file contents is shown.

A FITS file consists of one or more header-data units (HDUs),
concatenated one after the other. The `FITS` object therefore is
represented as a collection of these HDUs.

Iterate over HDUs in the file:
```jldoctest nasa_sample
julia> for hdu in f
           println(typeof(hdu))
       end
ImageHDU{Float32, 2}
ASCIITableHDU
```

### Image

Get information about the first image HDU:
```jldoctest nasa_sample
julia> hduimg = f[1]
File: docs/src/FOSy19g0309t_c2f.fits
HDU: 1
Mode: read-only
Type: Image
Datatype: Float32
Datasize: (2064, 2)
```

Each HDU can contain image data, or table data (either binary or
ASCII-formatted). For image extensions, get the size of the image
without reading it:
```jldoctest nasa_sample
julia> ndims(hduimg)
2

julia> size(hduimg)
(2064, 2)

julia> size(hduimg, 2)
2
```

Read an image from disk:
```jldoctest nasa_sample
julia> data = read(hduimg);  # read an image from disk

julia> data[1:4, :]
4×2 Matrix{Float32}:
 2.48511f-15  1.36115f-15
 1.56953f-15  1.10982f-15
 0.0          0.0
 1.12148f-15  9.71231f-16

julia> read(hduimg, 1:4, :) # read just a subset of image
4×2 Matrix{Float32}:
 2.48511f-15  1.36115f-15
 1.56953f-15  1.10982f-15
 0.0          0.0
 1.12148f-15  9.71231f-16
```

### Table

Show info about a binary table:
```jldoctest nasa_sample
julia> hdutable = f[2]
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
```

Read a column from the table:
```jldoctest nasa_sample
julia> read(hdutable, "CRVAL1")
2-element Vector{Float64}:
 1.0
 1.0

julia> read(hdutable, "CTYPE1")
2-element Vector{String}:
 "PIXEL"
 "PIXEL"
```

Table HDUs implement the [Tables.jl](https://tables.juliadata.org/stable/) interface, so you can load them into other table types, like [DataFrames](https://dataframes.juliadata.org/stable/).
```jldoctest nasa_sample
julia> using DataFrames

julia> df = DataFrame(hdutable);

julia> df[:, 1:5]
2×5 DataFrame
 Row │ CRVAL1   CRPIX1   CD1_1    DATAMIN  DATAMAX
     │ Float64  Float64  Float64  Float64  Float64
─────┼─────────────────────────────────────────────────
   1 │     1.0      1.0      1.0      0.0  2.73876e-15
   2 │     1.0      1.0      1.0      0.0  1.93483e-15
```

Variable length columns are not supported by the Tables.jl interface, and `Tables` methods will ignore them.

### Header

Read the entire header into memory and get values from it:
```jldoctest nasa_sample
julia> header = read_header(hduimg); # read the entire header from disk

julia> length(header)  # total number of records in header
162

julia> haskey(header, "NAXIS1")  # check if a key exists
true

julia> header["NAXIS1"]  # get value by keyword
2064

julia> header[4]  # get value by position
2064

julia> get_comment(header, "NAXIS")  # get comment for a given keyword
"Number of axes"
```


Read just a single header record without reading the entire header:
```jldoctest nasa_sample
julia> read_key(hduimg, 4)  # by position
("NAXIS1", 2064, "")

julia> read_key(hduimg, "NAXIS1")  # read by keyword
(2064, "")
```

Manipulate a header in memory:
```jldoctest nasa_sample
julia> header["NEWKEY"] = 10;  # change or add a keyword

julia> header["NEWKEY"]
10

julia> set_comment!(header, "NEWKEY", "this is a comment");

julia> get_comment(header, "NEWKEY")
"this is a comment"
```
Note that modifying the header does not update the file, unless this is explicitly written to the file.
```jldoctest nasa_sample
julia> haskey(read_header(hduimg), "NEWKEY")
false
```

### Close the file

```jldoctest nasa_sample
julia> close(f)
```
(`FITS` objects are also closed automatically when garbage collected.)

```@meta
DocTestFilters = nothing
```

## Writing to a file

Open a new file for writing:
```jldoctest example2
julia> using FITSIO

julia> f = FITS("newfile.fits", "w");
```
The second argument can be `"r"` (read-only; default), `"r+"`
(read-write) or `"w"` (write). In "write" mode, any existing file of
the same name is overwritten.

### Image

Write an image to the file:
```jldoctest example2
julia> data = reshape([1:10;], 2, 5)
2×5 Matrix{Int64}:
 1  3  5  7   9
 2  4  6  8  10

julia> write(f, data)  # Write a new image extension with the data

julia> close(f)
```
To write some header keywords in the new extension, pass a
`FITSHeader` instance as a keyword: `write(f, data; header=header)`


Overwrite image data in an existing file:
```jldoctest example2
julia> f = FITS("newfile.fits", "r+")  # Reopen the file in read-write mode
File: newfile.fits
Mode: "r+" (append)
HDUs: Num  Name  Type
      1          Image

julia> image_hdu = f[1];

julia> read(image_hdu)
2×5 Matrix{Int64}:
 1  3  5  7   9
 2  4  6  8  10
```

Write new data to the HDU

```jldoctest example2
julia> data = reshape([101:110;], 2, 5)
2×5 Matrix{Int64}:
 101  103  105  107  109
 102  104  106  108  110

julia> write(image_hdu, data)  # Overwrite the image

julia> read(image_hdu)
2×5 Matrix{Int64}:
 101  103  105  107  109
 102  104  106  108  110
```

### Table

Write a table to the file:
```jldoctest example2
julia> data = Dict("col1"=>[1., 2., 3.], "col2"=>[1, 2, 3])
Dict{String, Vector} with 2 entries:
  "col2" => [1, 2, 3]
  "col1" => [1.0, 2.0, 3.0]

julia> write(f, data)  # write a new binary table to a new extension

julia> hdutable = f[2]
File: newfile.fits
HDU: 2
Type: Table
Rows: 3
Columns: Name  Size  Type     TFORM
         col2        Int64    1K
         col1        Float64  1D

julia> read(hdutable, "col1")
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```

!!! tip "Compressed storage"
    Setting the file extension to `.gz` will automatically use GZIP compression and save on storage space.

    ```julia-repl
    julia> FITS("abc.fits", "w") do f # save the image uncompressed
               write(f, ones(200,200))
           end

    julia> filesize("abc.fits")
    325440

    julia> FITS("abc.fits.gz", "w") do f # save the image compressed
                write(f, ones(200,200))
           end

    julia> filesize("abc.fits.gz")
    2117
    ```

    Alternately the compression algorithm might be specified in square brackets after the filename. Check the [CFITSIO website](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/compression.html) for the details of this usage.

    ```julia-repl
    julia> FITS("abc.fits[compress R 100,100]", "w") do f # Rice algorithm with a 100 x 100 pixel tile size
               write(f, ones(200,200))
           end

    julia> filesize("abc.fits")
    8640
    ```

    !!! warn
        Compression is "loss-less" for images with integer pixel values, and might be lossy for floating-point images.
