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

`FITSIO.jl` can be installed using the built-in package manager:

```julia-repl
pkg> add FITSIO
```

```@meta
DocTestFilters = r"File: [a-zA-Z0-9/._]+"
```

## Quick start

The simplest way to write and read an image from a FITS file is by using the functions [`FITSIO.fitswrite`](@ref) and [`FITSIO.fitsread`](@ref).

```jldoctest
julia> fname, _ = mktemp(); # create a temporary file for this example

julia> using FITSIO

julia> FITSIO.fitswrite(fname, [1 2; 3 4]) # creates a new file with the data

julia> FITSIO.fitsread(fname) # opened in a read-only mode
2×2 Matrix{Int64}:
 1  2
 3  4

julia> FITSIO.fitsread(fname, 1, 1:2, 1) # read a section of HDU number 1
2-element Vector{Int64}:
 1
 3
```
This is not the most performant way, as the file is opened and closed on every invocation. However, this abstracts away the complexity.
!!! warn
    Currently, `fitswrite` overwrites an existing file, so this should be used with caution. In the future, this may allow
    appending to an existing file.

## Usage

In this example, we write to and read from a fits file.

We create a new temporary file using the `FITS` constructor:
```jldoctest example
julia> fname, _ = mktemp(); # create a temporary file for this example

julia> using FITSIO

julia> f = FITS(fname, "w") # create a new fits file
File: /tmp/jl_td8sL6
Mode: "w" (read-write)
No HDUs.

julia> length(f) # shows the number of HDUs in the file
0
```
In this example, we have used the file mode `"w"`, which will overwrite any existing file with the same name. To append to an existing file instead, we use the mode `"r+"`, and to open a file for reading, we may specify the mode `"r"`.

A FITS file consists of one or more header-data units (HDUs),
concatenated one after the other. The `FITS` object therefore is
represented as a collection of these HDUs. Each HDU can contain image data, or table data (either binary or
ASCII-formatted). Since we have just created a new file, there are no HDUs currently in the file.

### Image

We write an array to the file, which will automatically create an image HDU.

```jldoctest example
julia> data = Array(reshape(1:20, 4, 5))
4×5 Matrix{Int64}:
 1  5   9  13  17
 2  6  10  14  18
 3  7  11  15  19
 4  8  12  16  20

julia> write(f, data)

julia> length(f) # an HDU should now be added
1

julia> f # show information about the file
File: /tmp/jl_MdlNSQ
Mode: "w" (read-write)
HDUs: Num  Name  Type
      1          Image

julia> imghdu = f[1] # get the first image HDU
File: /tmp/jl_q95Lxk
HDU: 1
Mode: read-write
Type: Image
Datatype: Int64
Datasize: (4, 5)

julia> typeof(imghdu)
ImageHDU{Int64, 2}
```

We may query the properties of an image HDU analogous to an array.

```jldoctest example
julia> eltype(imghdu)
Int64

julia> ndims(imghdu)
2

julia> size(imghdu)
(4, 5)

julia> size(imghdu, 1)
4

julia> length(imghdu)
20
```

We may read the entire data array or a section of it using `read`. This is analogous to indexing into an array.

```jldoctest example
julia> read(imghdu)
4×5 Matrix{Int64}:
 1  5   9  13  17
 2  6  10  14  18
 3  7  11  15  19
 4  8  12  16  20

julia> read(imghdu, 1:3, 2:4)
3×3 Matrix{Int64}:
 5   9  13
 6  10  14
 7  11  15

julia> read(imghdu, 1:3, 3)
3-element Vector{Int64}:
  9
 10
 11

julia> read(imghdu, :, 2:3)
4×2 Matrix{Int64}:
 5   9
 6  10
 7  11
 8  12
```

### Header

We may read the header of the HDU as
```jldoctest example
julia> header = read_header(imghdu)
SIMPLE  =                    T / file does conform to FITS standard
BITPIX  =                   64 / number of bits per data pixel
NAXIS   =                    2 / number of data axes
NAXIS1  =                    4 / length of data axis 1
NAXIS2  =                    5 / length of data axis 2
EXTEND  =                    T / FITS dataset may contain extensions
COMMENT   FITS (Flexible Image Transport System) format is defined in 'Astronom
COMMENT   and Astrophysics', volume 376, page 359; bibcode: 2001A&A...376..359H

julia> typeof(header)
FITSHeader

julia> length(header) # number of keys
8
```
This reads the entire header from the HDU into memory. The header object behaves like a dictionary that may be queried using its keys:

```jldoctest example
julia> keys(header)
8-element Vector{String}:
 "SIMPLE"
 "BITPIX"
 "NAXIS"
 "NAXIS1"
 "NAXIS2"
 "EXTEND"
 "COMMENT"
 "COMMENT"

julia> haskey(header, "NAXIS")
true

julia> header["NAXIS"] # using the key name
2

julia> header[3] # using the key index
2
```

We may read the value and comment for a key as well, either using the key name or its index.

```jldoctest example
julia> read_key(imghdu, "NAXIS") # query by key name
(2, "number of data axes")

julia> read_key(imghdu, 3) # query by key index
("NAXIS", 2, "number of data axes")

julia> get_comment(header, "NAXIS")
"number of data axes"
```
We may modify the header in memory as
```jldoctest example
julia> header["NEWKEY"] = 10;  # change or add a keyword

julia> header["NEWKEY"]
10

julia> set_comment!(header, "NEWKEY", "this is a comment");

julia> get_comment(header, "NEWKEY")
"this is a comment"
```
!!! note
    Manipulating a header only changes it in memory until it is written to disk. The header object in memory is not connected to the fits file. To write some header keywords in the new extension, pass a [`FITSHeader`](@ref) instance as a keyword: `write(f, data; header=header)`

### Table

We append a binary table HDU to the file:
```jldoctest example
julia> data = Dict("col1"=>[1., 2., 3.], "col2"=>[1, 2, 3]);

julia> write(f, data)

julia> length(f) # check that a second HDU is now added
2

julia> f
File: /tmp/jl_q95Lxk
Mode: "w" (read-write)
HDUs: Num  Name  Type
      1          Image
      2          Table

julia> tablehdu = f[2]
File: /tmp/jl_q95Lxk
HDU: 2
Type: Table
Rows: 3
Columns: Name  Size  Type     TFORM
         col2        Int64    1K
         col1        Float64  1D

julia> typeof(tablehdu)
TableHDU
```

We may read a column from a table HDU as
```jldoctest example
julia> read(tablehdu, "col1")
3-element Vector{Float64}:
 1.0
 2.0
 3.0
```

Table HDUs implement the [Tables.jl](https://tables.juliadata.org/stable/) interface, so you can load them into other table types, like [DataFrames](https://dataframes.juliadata.org/stable/).

```jldoctest example
julia> using DataFrames

julia> df = DataFrame(tablehdu)
3×2 DataFrame
 Row │ col2   col1
     │ Int64  Float64
─────┼────────────────
   1 │     1      1.0
   2 │     2      2.0
   3 │     3      3.0
```
Variable length columns are not supported by the `Tables.jl` interface, and `Table`s methods will ignore them.

### Iterating over the HDUs

We may iterate over the fits file to access individual HDUs.
```jldoctest example
julia> for hdu in f
           println(typeof(hdu))
       end
ImageHDU{Int64, 2}
TableHDU
```

### Write to disk

Close the file to write the in-memory FITS file to disk.

```jldoctest example
julia> close(f)
```
FITS objects are also closed automatically when garbage collected.

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
