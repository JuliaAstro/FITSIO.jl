module FITSIO

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
       copy_section,
       default_header

import Base: getindex,
             setindex!,
             length,
             show,
             read,
             read!,
             write,
             close,
             ndims,
             size,
             get,
             getkey,
             haskey,
             keys,
             values,
             eltype,
             deleteat!,
             delete!

# Deal with compatibility issues.
using Printf
import Base: iterate, lastindex
# Libcfitsio submodule

## Have to manually import while deprecating `libcfitsio_version`.
## After removing that, can return to using CFITSIO
import CFITSIO: FITSFile,
                FITSMemoryHandle,
                fits_open_file,
                fits_open_diskfile,
                fits_create_file,
                fits_create_diskfile,
                fits_assert_open,
                fits_file_mode,
                fits_create_img,
                fits_close_file,
                fits_write_pix,
                fits_get_num_hdus,
                fits_movabs_hdu,
                fits_delete_hdu,
                fits_get_img_size,
                fits_get_img_type,
                fits_get_img_equivtype,
                type_from_bitpix,
                fits_read_pix,
                fits_read_subset,
                fits_copy_image_section,
                fits_file_name,
                fits_get_img_dim,
                fits_read_tdim,
                fits_write_tdim,
                fits_read_col,
                fits_write_col,
                fits_get_num_cols,
                fits_get_num_rows,
                fits_get_colnum,
                fits_get_eqcoltype,
                fits_read_descript,
                fits_update_key,
                fits_write_comment,
                fits_write_history,
                fits_get_hdrspace,
                fits_hdr2str,
                fits_read_keyword,
                fits_read_keyn,
                fits_open_memfile,
                fits_movnam_hdu,
                fits_get_hdu_num,
                fits_get_hdu_type

# There are a few direct `ccall`s to libcfitsio in this module. For this, we
# need a few non-exported things from Libcfitsio: the shared library handle,
# and a helper function for raising errors.
import CFITSIO: libcfitsio,
                fits_assert_ok,
                fits_assert_isascii

import CFITSIO
@deprecate libcfitsio_version CFITSIO.libcfitsio_version

import Tables

## DEPRECATED
module Libcfitsio
    using Reexport
    @reexport using CFITSIO
end

# HDU Types
abstract type HDU end

mutable struct ImageHDU{T<:Real,N} <: HDU
    fitsfile::FITSFile
    ext::Int

    function ImageHDU(fitsfile::FITSFile, ext::Integer)
        fits_assert_open(fitsfile)
        fits_movabs_hdu(fitsfile, ext)
        N = Int(fits_get_img_dim(fitsfile))
        bitpix = fits_get_img_equivtype(fitsfile)
        T = type_from_bitpix(bitpix)
        new{T,N}(fitsfile, ext)
    end
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
    FITS(filename::String[, mode::String = "r"]; extendedparser = true)

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

The keyword argument `extendedparser` may be used to enable or disable the
[extended filename parser](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html).
If disabled, `filename` is treated exactly as the name of the file and is not tokenized into
parameters.
"""
mutable struct FITS
    fitsfile::FITSFile
    filename::AbstractString
    mode::AbstractString
    hdus::Dict{Int, HDU}

    # Hold on to memory if backed by a julia buffer
    memhandle::FITSMemoryHandle
    hidden::Any

    function FITS(filename::AbstractString, mode::AbstractString="r"; extendedparser = true)
        openfn = extendedparser ? fits_open_file : fits_open_diskfile
        createfn = extendedparser ? fits_create_file : fits_create_diskfile
        f = (mode == "r"                      ? openfn(filename, 0)    :
             mode == "r+" && isfile(filename) ? openfn(filename, 1)    :
             mode == "r+"                     ? createfn(filename)     :
             mode == "w"                      ?
                (rm(filename, force = true); createfn(filename)) :
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
    FITS(f::Function, args...; kwargs...)

Apply the function `f` to the result of `FITS(args...; kwargs...)` and close the
resulting file descriptor upon completion.
"""
function FITS(f::Function, args...; kwargs...)
    io = FITS(args...; kwargs...)
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

end # module
