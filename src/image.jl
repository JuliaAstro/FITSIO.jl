# ImageHDU methods

# Display the image datatype and dimensions
function show(io::IO, hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    if isdeleted(hdu)
        print(io, "Deleted HDU")
        return
    end
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_type(hdu.fitsfile)
    equivbitpix = fits_get_img_equivtype(hdu.fitsfile)
    sz = fits_get_img_size(hdu.fitsfile)
    mode = fits_file_mode(hdu.fitsfile)

    if bitpix == equivbitpix
        datainfo = string(type_from_bitpix(equivbitpix))
    else
        datainfo = @sprintf("%s (physical: %s)",
                            type_from_bitpix(equivbitpix),
                            type_from_bitpix(bitpix))
    end

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(hdu.ext)$(fits_get_ext_info_string(hdu.fitsfile))
    Mode: $(mode == 0 ? "read-only" : "read-write")
    Type: Image
    Datatype: $datainfo
    Datasize: $(Tuple(sz))""")
end

# Get image dimensions
"""
    ndims(hdu::ImageHDU)

Get number of image dimensions, without reading the image into memory.
"""
ndims(::ImageHDU{<:Any,N}) where N = N

"""
    size(hdu::ImageHDU)
    size(hdu::ImageHDU, i)

Get image dimensions (or `i`th dimension), without reading the image
into memory.
"""
function size(hdu::ImageHDU{<:Any,N}) where N
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile, Val(N))
    NTuple{N,Int}(sz)
end

size(hdu::ImageHDU, i::Integer) = size(hdu)[i]

"""
    length(hdu::ImageHDU)

Get total number of pixels in image (product of `size(hdu)`).
"""
length(hdu::ImageHDU) = prod(size(hdu))

# `lastindex` is needed so that hdu[:] can throw DimensionMismatch
# when ndim != 1, rather than no method.
lastindex(hdu::ImageHDU) = length(hdu)

"""
    eltype(hdu::ImageHDU)

Return the element type of the image in `hdu`.
"""
eltype(::ImageHDU{T}) where T = T

"""
    fitsread(filename::AbstractString[, hduindex = 1[, arrayindices...]]; extendedparser = true)

Convenience function to read in an image corresponding to the HDU at index `hduindex` contained in
the FITS file named `filename`.
If `arrayindices` are provided, only a slice of the image corresponding to the indices is read in.

Functionally `fitsread(filename, hduindex, arrayindices...; extendedparser)` is equivalent to

```julia
FITS(filename, "r"; extendedparser = extendedparser) do f
    read(f[hduindex], arrayindices...)
end
```

The keyword argument `extendedparser` may be used to enable or disable the
[extended filename parser](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html).
If disabled, `filename` is treated exactly as the name of the file and is not tokenized into
parameters.

!!! note
    Julia follows a column-major array indexing convention, so the indices provided must account for this.
    In particular this means that FITS files created externally following a row-major convention (eg. using astropy)
    will have the sequence of axes flipped when read in using FITSIO.

See also: [`read`](@ref)
"""
function fitsread(filename::AbstractString, hduindex = 1, arrayindices...; extendedparser = true)
    FITS(filename, "r"; extendedparser = extendedparser) do f
        read(f[hduindex], arrayindices...)
    end
end

# Read a full image from an HDU
"""
    read(hdu::ImageHDU)
    read(hdu::ImageHDU, range...)

Read the data array or a subset thereof from disk. The first form
reads the entire data array. The second form reads a slice of the array
given by the specified ranges or integers. Dimensions specified by integers will be
dropped in the returned array, while those specified by ranges will be retained.

!!! note
    Julia follows a column-major array indexing convention, so the indices provided must account for this.
    In particular this means that FITS files created externally following a row-major convention (eg. using astropy)
    will have the sequence of axes flipped when read in using FITSIO.
"""
function read(hdu::ImageHDU)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    data = Array{eltype(hdu)}(undef, size(hdu))
    fits_read_pix(hdu.fitsfile, data)
    data
end
read(hdu::ImageHDU{<:Any,0}) = (assert_open(hdu); nothing)

#= Inplace read
This requires a contiguous array. Lacking a type to dispatch upon,
we rely on runtime checks.
=#

iscontiguous(array::Union{Array,StridedArray{<:Any,0}}) = true
function iscontiguous(array)
    strd = strides(array)
    sz = size(array)
    isone(strd[1]) && checkcontiguous(Base.tail(strd),sz[1],Base.tail(sz)...)
end

function checkcontiguous(strd,s,d,sz...)
    strd[1] == s && checkcontiguous(Base.tail(strd),s*d,sz...)
end
checkcontiguous(::Tuple{},args...) = true

"""
    read!(hdu::ImageHDU, A::StridedArray)
    read!(hdu::ImageHDU, A::StridedArray, range...)

Read the data or a subset thereof from disk, and save it in a
pre-allocated output array `A`.
The first form reads the entire data from disk.
The second form reads a slice of the array given by the specified ranges or integers.
The array `A` needs to have the same length as the number of elements to be read in.
Additionally `A` needs to be stored contiguously in memory.

!!! note
    Julia follows a column-major array indexing convention, so the indices provided must account for this.
    In particular this means that FITS files created externally following a row-major convention (eg. using astropy)
    will have the sequence of the axes flipped when read in using FITSIO.
"""
function read!(hdu::ImageHDU, array::StridedArray{<:Real})
    read!(hdu, reshape(array, size(hdu)))
end
function read!(hdu::ImageHDU{<:Real,N}, array::StridedArray{<:Real,N}) where {N}

    if !iscontiguous(array)
        throw(ArgumentError("the output array needs to be contiguous"))
    end

    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    if length(hdu) != length(array)
        throw(DimensionMismatch("length of the output array does not match the number of elements to be read in"))
    end

    fits_read_pix(hdu.fitsfile, array)
    return array
end

# _checkbounds methods copied from Julia v0.4 Base.
_checkbounds(sz, i::Integer) = 1 <= i <= sz
_checkbounds(sz, i::Colon) = true
_checkbounds(sz, r::AbstractRange{<:Integer}) =
    (isempty(r) || (minimum(r) >= 1 && maximum(r) <= sz))

# helper functions for constructing cfitsio indexing vectors in read(hdu, ...)
_first(i::Union{Integer, AbstractRange}) = first(i)
_first(::Colon) = 1
_last(sz, i::Union{Integer, AbstractRange}) = last(i)
_last(sz, ::Colon) = sz
_step(r::AbstractRange) = step(r)
_step(::Union{Integer, Colon}) = 1

# Shape of array to create for read(hdu, ...), dropping trailing
# scalars. This is simpler than in Base because we are guaranteed that
# length(I) == length(sz).
@inline _index_shape(sz, I...) = _index_shape_dim(sz, 1, I...)
_index_shape_dim(sz, dim, I::Integer...) = ()
_index_shape_dim(sz, dim, ::Colon) = (sz[dim],)
_index_shape_dim(sz, dim, r::AbstractRange) = (length(r),)
@inline _index_shape_dim(sz, dim, ::Colon, I...) =
    tuple(sz[dim], _index_shape_dim(sz, dim+1, I...)...)
@inline _index_shape_dim(sz, dim, ::Integer, I...) =
    _index_shape_dim(sz, dim+1, I...)
@inline _index_shape_dim(sz, dim, r::AbstractRange, I...) =
    tuple(length(r), _index_shape_dim(sz, dim+1, I...)...)

# Read a subset of an ImageHDU
function read_internal(hdu::ImageHDU, I::Union{AbstractRange{<:Integer}, Integer, Colon}...)

    # check number of indices and bounds. Note that number of indices and
    # array dimension must match, unlike in Arrays. Array-like behavior could
    # be supported in the future with care taken in constructing first, last,
    if length(I) != ndims(hdu)
        throw(DimensionMismatch("number of indices must match dimensions"))
    end

    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    sz = size(hdu)
    for i=1:ndims(hdu)
        _checkbounds(sz[i], I[i]) || throw(BoundsError())
    end

    # construct first, last and step vectors
    firsts = Clong[_first(idx) for idx in I]
    lasts = Clong[_last(sz[i], I[i]) for i=1:length(sz)]
    steps = Clong[_step(idx) for idx in I]

    # construct output array
    data = Array{eltype(hdu)}(undef, _index_shape(sz, I...))

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

function read_internal!(hdu::ImageHDU, array::StridedArray,
    I::Union{AbstractRange{<:Integer}, Integer, Colon}...)

    if !iscontiguous(array)
        throw(ArgumentError("the output array needs to be contiguous"))
    end
    # check number of indices and bounds. Note that number of indices and
    # array dimension must match, unlike in Arrays. Array-like behavior could
    # be supported in the future with care taken in constructing first, last,
    if length(I) != ndims(hdu)
        throw(DimensionMismatch("number of indices must match dimensions"))
    end

    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    sz = size(hdu)
    for i = 1:ndims(hdu)
        _checkbounds(sz[i], I[i]) || throw(BoundsError())
    end

    ninds = _index_shape(sz, I...)

    if length(array) != prod(ninds)
        throw(DimensionMismatch("length of output array does not match the number of elements to be read in"))
    end

    # construct first, last and step vectors
    firsts = Clong[_first(idx) for idx in I]
    lasts = Clong[_last(sz[i], I[i]) for i=1:ndims(hdu)]
    steps = Clong[_step(idx) for idx in I]

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, array)
    array
end

# general method and version that returns a single value rather than 0-d array
read(hdu::ImageHDU, I::Union{AbstractRange{<:Integer}, Integer, Colon}...) =
    read_internal(hdu, I...)
read(hdu::ImageHDU, I::Integer...) = read_internal(hdu, I...)[1]

read!(hdu::ImageHDU, array::StridedArray{<:Real}, I::Union{AbstractRange{<:Integer}, Integer, Colon}...) =
    read_internal!(hdu, array, I...)
read!(hdu::ImageHDU, array::StridedArray{<:Real}, I::Integer...) = read_internal!(hdu, array, I...)[1]

"""
    fitswrite(filename::AbstractString, data; extendedparser = true, kwargs...)

Convenience function to write the image array `data` to a file named `filename`.

Functionally `fitswrite(filename, data; extendedparser, kwargs...)` is equivalent to

```julia
FITS(filename, "w"; extendedparser = extendedparser) do f
    write(f, data; kwargs...)
end
```

The keyword argument `extendedparser` may be used to enable or disable the
[extended filename parser](https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node83.html).
If disabled, `filename` is treated exactly as the name of the file and is not tokenized into
parameters.

!!! warn "Warning"
    Existing files with the same name will be overwritten.

See also: [`write`](@ref)
"""
function fitswrite(filename::AbstractString, data; extendedparser = true, kwargs...)
    FITS(filename, "w", extendedparser = extendedparser) do f
        write(f, data; kwargs...)
    end
end

"""
    write(f::FITS, data::StridedArray{<:Real};
        header=nothing, name=nothing, ver=nothing, checksum::Bool=false)

Add a new image HDU to FITS file `f` with contents `data`. The
following array element types are supported: `UInt8`, `Int8`,
`UInt16`, `Int16`, `UInt32`, `Int32`, `Int64`, `Float32`,
`Float64`. If a `FITSHeader` object is passed as the `header` keyword
argument, the header will also be added to the new HDU.

The keyword argument `checksum` may be set to `true` to enable
writing a checksum for the HDU. This is not enabled by default, as it
can significantly slow down writing large HDUs.

!!! note "Contiguous data"
    The `data` array must be stored contiguously in memory.
"""
function write(f::FITS, data::StridedArray{<:Real};
               header::Union{Nothing, FITSHeader}=nothing,
               name::Union{Nothing, String}=nothing,
               ver::Union{Nothing, Integer}=nothing,
               checksum::Bool=false)

    if !iscontiguous(data)
        throw(ArgumentError("data to be written out needs to be contiguous"))
    end

    if f.mode == "r"
        throw(ArgumentError("FITS file has been opened in read-only mode"))
    end

    ndims(data) == 0 && throw(ArgumentError("data must have at least one dimension"))

    fits_create_img(f.fitsfile, data)

    if isa(header, FITSHeader)
        fits_write_header(f.fitsfile, header, true)
    end
    if isa(name, String)
        fits_update_key(f.fitsfile, "EXTNAME", name)
    end
    if isa(ver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", ver)
    end
    fits_write_pix(f.fitsfile, data)
    checksum && fits_write_chksum(f.fitsfile)
    nothing
end

"""
    write(hdu::ImageHDU, data::StridedArray{<:Real}; checksum::Bool=false)

Write data to an existing image HDU.

The keyword argument `checksum` may be set to `true` to enable
writing a checksum for the HDU. This is not enabled by default, as it
can significantly slow down writing large HDUs.

!!! note
    The `data` array must be stored contiguously in memory.
"""
function write(hdu::ImageHDU, data::StridedArray{<:Real}; checksum::Bool=false)

    if !iscontiguous(data)
        throw(ArgumentError("data to be written out needs to be contiguous"))
    end

    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    if fits_file_mode(hdu.fitsfile) == 0
        throw(ArgumentError("FITS file has been opened in read-only mode"))
    end

    # Ensure sizes are equal
    hdu_size = size(hdu)
    data_size = size(data)

    if hdu_size != data_size
        error("size of HDU $(hdu_size) not equal to size of data $(data_size).")
    end

    fits_write_pix(hdu.fitsfile, data)
    checksum && fits_write_chksum(hdu.fitsfile)
    nothing
end

# Copy a rectangular section of an image and write it to a new FITS
# primary image or image extension. The new image HDU is appended to
# the end of the output file; all the keywords in the input image will
# be copied to the output image. The common WCS keywords will be
# updated if necessary to correspond to the coordinates of the section.

# convert a range to a string that cfitsio understands
cfitsio_range_string(r::UnitRange) = @sprintf "%d:%d" first(r) last(r)
cfitsio_range_string(r::StepRange) =
    @sprintf "%d:%d:%d" first(r) last(r) step(r)

"""
    copy_section(hdu, dest, r...)

Copy a rectangular section of an image and write it to a new FITS
primary image or image extension in `FITS` object `dest`. The new
image HDU is appended to the end of `dest`. All the keywords in the
input image will be copied to the output image. The common WCS
keywords will be updated if necessary to correspond to the coordinates
of the section.

# Examples

Copy the lower-left 200 x 200 pixel section of the image in `hdu`
to an open file, `f`

```
copy_section(hdu, f, 1:200, 1:200)
```

Same as above but only copy odd columns in y:

```
copy_section(hdu, f, 1:200, 1:2:200)
```
"""
function copy_section(hdu::ImageHDU, dest::FITS, r::AbstractRange{Int}...)
    assert_open(hdu)
    fits_assert_open(dest.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_copy_image_section(hdu.fitsfile, dest.fitsfile,
                            join([cfitsio_range_string(ri) for ri in r], ','))
end
