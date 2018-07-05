__precompile__()


module FITSIO

using Printf

export FITS,
       HDU,
       ImageHDU,
       TableHDU,
       ASCIITableHDU,
       FITSHeader,
       read_key,
       write_key,
       read_header,
       get_comment,
       set_comment!,
       copy_section

import Base: getindex,
             setindex!,
             length,
             show,
             read,
             write,
             close,
             ndims,
             size,
             endof,
             haskey,
             keys,
             values,
             start,
             next,
             done

# Libcfitsio submodule
include("libcfitsio.jl")

using .Libcfitsio

# There are a few direct `ccall`s to libcfitsio in this module. For this, we
# need a few non-exported things from Libcfitsio: the shared library handle,
# and a helper function for raising errors. TYPE_FROM_BITPIX is awkwardly
# defined in Libcfitsio, even though it is not used there.
import .Libcfitsio: libcfitsio,
                    fits_assert_ok,
                    fits_assert_isascii,
                    TYPE_FROM_BITPIX

# HDU Types
abstract type HDU end

mutable struct ImageHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

mutable struct TableHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

mutable struct ASCIITableHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

# FITS
#
# The FITS type represents a FITS file. It holds a reference to a
# FITSFile object (basically the low-level CFITSIO pointer). It also
# holds a reference to all of the previously accessed HDU
# objects. This is so that only a single HDU object is created for
# each extension in the file. It also allows a FITS object to tell
# previously created HDUs about events that happen to the file, such
# as deleting extensions. This could be done by, e.g., setting ext=-1
# in the HDU object.
"""
    FITS(filename::String, mode::String="r")

Open or create a FITS file. `mode` can be one of `"r"` (read-only),
`"r+"` (read-write) or `"w"` (write). In "write" mode, any existing
file of the same name is overwritten.

A `FITS` object is a collection of "Header-Data Units" (HDUs) and
supports the following operations:

* `f[i]`: Return the `i`-th HDU.

* `f[name]` or `f[name, ver]`: Return the HDU containing the given the
  given EXTNAME (or HDUNAME) keyword (a String), and optionally the
  given EXTVER (or HDUVER) number (an Integer).

* Iteration:
  ```
  for hdu in f
      ...
  end
  ```
"""
FITS
mutable struct FITS
    fitsfile::FITSFile
    filename::AbstractString
    mode::AbstractString
    hdus::Dict{Int, HDU}

    # Hold on to memory if backed by a julia buffer
    memhandle::FITSMemoryHandle
    hidden::Any

    function FITS(filename::AbstractString, mode::AbstractString="r")
        f = (mode == "r"                      ? fits_open_file(filename, 0)    :
             mode == "r+" && isfile(filename) ? fits_open_file(filename, 1)    :
             mode == "r+"                     ? fits_create_file(filename)     :
             mode == "w"                      ? fits_create_file("!"*filename) :
             error("invalid open mode: $mode"))

        new(f, filename, mode, Dict{Int, HDU}(), FITSMemoryHandle(), nothing)
    end

    function FITS(data::Vector{UInt8}, mode::AbstractString="r", filename = "")
        @assert mode == "r"
        f, handle = fits_open_memfile(data, 0)
        new(f, filename, mode, Dict{Int, HDU}(), handle, data)
    end
end

"""
    FITS(f::Function, args...)

Apply the function `f` to the result of `FITS(args...)` and close the
resulting file descriptor upon completion.
"""
function FITS(f::Function, args...)
    io = FITS(args...)
    try
        f(io)
    finally
        close(io)
    end
end


"""
    FITSHeader(keys::Vector{String}, values::Vector, comments::Vector{String})

An in-memory representation of the header of an HDU. It stores the
(key, value, comment) information for each 80-character "card" in a header.

Note that this structure is not linked to a FITS file in any way; it
is just a convenient structure for storing the header contents after
reading from a file. (This is similar to how an `Array` returned by
`read(f[1])` is not linked to the FITS file `f`.)  Manipulating a
`FITSHeader` will therefore have no immediate impact on any file, even
if it was created by `read_header(::HDU)`.  You can, however, write a
`FITSHeader` to a file using the `write(::FITS, ...)` methods that
append a new HDU to a file.
"""
mutable struct FITSHeader
    keys::Vector{String}
    values::Vector{Any}
    comments::Vector{String}
    map::Dict{String, Int}

    function FITSHeader(keys::Vector{String}, values::Vector,
                        comments::Vector{String})
        if ((length(keys) != length(values)) ||
            (length(keys) != length(comments)))
            error("keys, values, comments must be same length")
        end
        map = Dict{String, Int}()
        for i in 1:length(keys)
          map[keys[i]] = i
        end
        new(keys, Vector{Any}(values), comments, map)
    end
end

include("fits.jl")  # FITS methods
include("header.jl")  # FITSHeader methods
include("image.jl")  # ImageHDU methods
include("table.jl")  # TableHDU & ASCIITableHDU methods

"""
    libcfitsio_version() -> VersionNumber

Return the version of the underlying CFITSIO library

# Example

```julia
julia> FITSIO.libcfitsio_version()
v"3.37.0"
```

"""
function libcfitsio_version(version=fits_get_version())
    # fits_get_version returns a float. e.g., 3.341f0. We parse that
    # into a proper version number. E.g., 3.341 -> v"3.34.1"
    v = Int(round(1000 * version))
    x = div(v, 1000)
    y = div(rem(v, 1000), 10)
    z = rem(v, 10)
    VersionNumber(x, y, z)
end

end # module
