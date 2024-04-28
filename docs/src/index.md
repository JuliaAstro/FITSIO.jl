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

To open an existing file for reading:
```julia-repl
julia> using FITSIO

julia> f = FITS("file.fits")
File: file.fits
Mode: "w" (read-write)
HDUs: Num  Name  Type   
      1          Image  
      2          Table  
```
(At the REPL, information about the file contents is shown.)


A FITS file consists of one or more header-data units (HDUs),
concatenated one after the other. The `FITS` object therefore is
represented as a collection of these HDUs.

Get information about the first HDU:
```julia-repl
julia> f[1]
File: file.fits
HDU: 1
Type: Image
Datatype: Float64
Datasize: (800, 800)
```

Iterate over HDUs in the file:
```julia-repl
julia> for hdu in f; println(typeof(hdu)); end
FITSIO.ImageHDU
FITSIO.TableHDU
```


Each HDU can contain image data, or table data (either binary or
ASCII-formatted). For image extensions, get the size of the image
without reading it:
```julia-repl
julia> ndims(f[1])
    2

julia> size(f[1])
(800,800)

julia> size(f[1], 2)
800
```


Read an image from disk:
```julia-repl
julia> data = read(f[1]);  # read an image from disk

julia> data = read(f[1], :, 790:end);  # read just a subset of image
```


Show info about a binary table:
```julia-repl
julia> f[2]
File: file.fits
HDU: 2
Type: Table
Rows: 20
Columns: Name  Size  Type    TFORM  
         col2        String  5A     
         col1        Int64   1K     
```


Read a column from the table:
```julia-repl
julia> data = read(f[2], "col1")
```

Table HDUs implement the [Tables.jl](https://tables.juliadata.org/stable/) interface, so you can load them into other table types, like [DataFrames](https://dataframes.juliadata.org/stable/).
```julia-repl
julia> df = DataFrame(f[2])
```
Variable length columns are not supported by the Tables.jl interface, and `Tables` methods will ignore them.


Read the entire header into memory and get values from it:
```julia-repl
julia> header = read_header(f[1]);  # read the entire header from disk

julia> length(header)  # total number of records in header
17

julia> haskey(header, "NAXIS1")  # check if a key exists
true

julia> header["NAXIS1"]  # get value by keyword
800

julia> header[4]  # get value by position
800

julia> get_comment(header, "NAXIS")  # get comment for a given keyword
"length of data axis 1"
```


Read just a single header record without reading the entire header:
```julia-repl
julia> read_key(f[1], 4)  # by position
("NAXIS1",800,"length of data axis 1")

julia> read_key(f[1], "NAXIS1")  # read by keyword
(800,"length of data axis 1")
```


Manipulate a header in memory:
```julia-repl
julia> header["NEWKEY"] = 10  # change or add a keyword

julia> set_comment!(header, "NEWKEY", "this is a comment")
```


Close the file:
```julia-repl
julia> close(f)
```
(`FITS` objects are also closed automatically when garbage collected.)


Open a new file for writing:
```julia-repl
julia> f = FITS("newfile.fits", "w");
```
The second argument can be `"r"` (read-only; default), `"r+"`
(read-write) or `"w"` (write). In "write" mode, any existing file of
the same name is overwritten.


Write an image to the file:
```julia-repl
julia> data = reshape([1:100;], 5, 20)

julia> write(f, data)  # Write a new image extension with the data
julia> close(f)
```
To write some header keywords in the new extension, pass a
`FITSHeader` instance as a keyword: `write(f, data; header=header)`


Overwrite image data in an existing file:
```julia-repl
julia> f = FITS("newfile.fits", "r+")  # Reopen the file in read-write mode
julia> data = reshape([101:200;], 5, 20)  # Prepare new image data
julia> image_hdu = f[1]
julia> write(image_hdu, data)  # Overwrite the image
```


Write a table to the file:
```julia-repl
julia> data = Dict("col1"=>[1., 2., 3.], "col2"=>[1, 2, 3]);

julia> write(f, data)  # write a new binary table to a new extension
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
