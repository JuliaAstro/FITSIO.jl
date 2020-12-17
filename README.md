FITSIO.jl
=========

Flexible Image Transport System (FITS) support for Julia


[![Build Status](https://github.com/juliaastro/FITSIO.jl/workflows/CI/badge.svg?branch=master)](https://github.com/juliaastro/FITSIO.jl/actions)
[![Coverage Status](http://img.shields.io/coveralls/JuliaAstro/FITSIO.jl.svg?style=flat-square)](https://coveralls.io/r/JuliaAstro/FITSIO.jl?branch=master)


[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](http://juliaastro.github.io/FITSIO.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliahci.github.io/FITSIO.jl/dev)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Note:** The `Libcfitsio` submodule has been moved to [CFITSIO.jl](https://github.com/JuliaAstro/CFITSIO.jl) and will be deprecated in a future release.


## Usage

For more in-depth usage and examples, see [the documentation](http://juliaastro.github.io/FITSIO.jl/stable/)

```julia

julia> using FITSIO

julia> f = FITS("file.fits")
File: file.fits
Mode: "w" (read-write)
HDUs: Num  Name  Type
      1          Image
      2          Table

julia> f[1]
File: file.fits
HDU: 1
Type: Image
Datatype: Float64
Datasize: (800, 800)

# read an image from disk
julia> data = read(f[1]);  

# read just a subset of image
julia> data = read(f[1], :, 790:end);  

# Get info about binary table
julia> f[2]
File: file.fits
HDU: 2
Type: Table
Rows: 20
Columns: Name  Size  Type    TFORM  
         col2        String  5A     
         col1        Int64   1K     

# Read a column from the table:
 julia> data = read(f[2], "col1")

# Read the entire header into memory
julia> header = read_header(f[1]);

julia> header["NAXIS1"]  # get value by keyword
800

julia> header[4]  # get value by position
800

# Read single keys into memory
julia> read_key(f[1], 4) # by position
("NAXIS1",800,"length of data axis 1")

julia> read_key(f[1], "NAXIS1") # by keyword
(800,"length of data axis 1")

# write data to file
julia> FITS("new_file.fits", "w") do f
           write(f, data)
       end
```
