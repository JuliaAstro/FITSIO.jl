# ImageHDU methods

# Display the image datatype and dimensions
function show(io::IO, hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_type(hdu.fitsfile)
    equivbitpix = fits_get_img_equivtype(hdu.fitsfile)
    sz = fits_get_img_size(hdu.fitsfile)

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
    Type: Image
    Datatype: $datainfo
    Datasize: $(tuple(sz...))""")
end

# Get image dimensions
"""
    ndims(hdu::ImageHDU)

Get number of image dimensions, without reading the image into memory.
"""
function ndims(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_get_img_dim(hdu.fitsfile)
end

"""
    size(hdu::ImageHDU)
    size(hdu::ImageHDU, i)

Get image dimensions (or `i`th dimension), without reading the image
into memory.
"""
function size(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    tuple(sz...)
end

function size(hdu::ImageHDU, i::Integer)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    sz[i]
end

"""
    length(hdu::ImageHDU)

Get total number of pixels in image (product of ``size(hdu)``).
"""
length(hdu::ImageHDU) = prod(size(hdu))

# `lastindex` is needed so that hdu[:] can throw DimensionMismatch
# when ndim != 1, rather than no method.
lastindex(hdu::ImageHDU) = length(hdu::ImageHDU)

"""
    eltype(hdu::ImageHDU)

Return the element type of the image in `hdu`.
"""
function eltype(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    type_from_bitpix(bitpix)
end

# Read a full image from an HDU
"""
    read(hdu::ImageHDU)
    read(hdu::ImageHDU, range...)

Read the data array or a subset thereof from disk. The first form
reads the entire data array. The second form reads a slice of the array
given by the specified ranges or integers.
"""
function read(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    data = Array{type_from_bitpix(bitpix)}(undef, sz...)
    fits_read_pix(hdu.fitsfile, data)
    data
end

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
    read!(hdu::ImageHDU, array::StridedArray)
    read!(hdu::ImageHDU, array::StridedArray, range...)

Read the data or a subset thereof from disk, and save it in a 
pre-allocated output array.
The first form reads the entire data from disk. 
The second form reads a slice of the array given by the specified ranges or integers.
The output array needs to have the same shape as the data range to be read in. 
Additionally the output array needs to be contiguously stored in memory.
"""
function read!(hdu::ImageHDU, array::StridedArray{T,N}) where {T,N}

    if !iscontiguous(array)
        throw(ArgumentError("output array needs to be contiguously stored"))
    end

    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    
    if type_from_bitpix(bitpix) != T
        throw(TypeError(:read!,"",type_from_bitpix(bitpix),T))
    end
    
    if ndims(hdu) != N
        throw(DimensionMismatch("array dimensions do not match the data. "*
        "Data has $(ndims(hdu)) dimensions whereas the output array has $N dimensions."))
    end
    
    # Maybe this can be a ShapeMismatch when there's a decision on #16717
    if Tuple(sz) != size(array)
        throw(DimensionMismatch("size of the array does not "*
        "match the data. Data has a size of $(Tuple(sz)) whereas the output array "*
        "has a size of $(size(array))"))
    end

    fits_read_pix(hdu.fitsfile, array)
    array
end

# _checkbounds methods copied from Julia v0.4 Base.
_checkbounds(sz, i::Integer) = 1 <= i <= sz
_checkbounds(sz, i::Colon) = true
_checkbounds(sz, r::AbstractRange{Int}) =
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
function read_internal(hdu::ImageHDU, I::Union{AbstractRange{Int}, Integer, Colon}...)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)

    # check number of indices and bounds. Note that number of indices and
    # array dimension must match, unlike in Arrays. Array-like behavior could
    # be supported in the future with care taken in constructing first, last,
    if length(I) != length(sz)
        throw(DimensionMismatch("number of indices must match dimensions"))
    end
    for i=1:length(sz)
        _checkbounds(sz[i], I[i]) || throw(BoundsError())
    end

    # construct first, last and step vectors
    firsts = Clong[_first(idx) for idx in I]
    lasts = Clong[_last(sz[i], I[i]) for i=1:length(sz)]
    steps = Clong[_step(idx) for idx in I]

    # construct output array
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    data = Array{type_from_bitpix(bitpix)}(undef, _index_shape(sz, I...))

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

function read_internal!(hdu::ImageHDU, array::StridedArray{T,N}, 
    I::Union{AbstractRange{Int}, Integer, Colon}...) where {T,N}

    if !iscontiguous(array)
        throw(ArgumentError("output array needs to be contiguously stored"))
    end

    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # check that the output array has the right type
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    if type_from_bitpix(bitpix) != T
        throw(TypeError(:read!,"",type_from_bitpix(bitpix),T))
    end
    
    sz = fits_get_img_size(hdu.fitsfile)

    # check number of indices and bounds. Note that number of indices and
    # array dimension must match, unlike in Arrays. Array-like behavior could
    # be supported in the future with care taken in constructing first, last,
    if length(I) != length(sz)
        throw(DimensionMismatch("number of indices must match dimensions"))
    end

    for i=1:length(sz)
        _checkbounds(sz[i], I[i]) || throw(BoundsError())
    end

    ninds = _index_shape(sz, I...)
    if length(ninds) != N
        throw(DimensionMismatch("number of dimensions to be read must match the "*
            "number of dimensions of the output array"))
    end

    if size(array) != ninds
        throw(DimensionMismatch("size of the data slice must match that of the output array. "*
            "Data has a size of $ninds whereas the output array has a size of $(size(output))"))
    end

    # construct first, last and step vectors
    firsts = Clong[_first(idx) for idx in I]
    lasts = Clong[_last(sz[i], I[i]) for i=1:length(sz)]
    steps = Clong[_step(idx) for idx in I]

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, array)
    array
end

# general method and version that returns a single value rather than 0-d array
read(hdu::ImageHDU, I::Union{AbstractRange{Int}, Int, Colon}...) =
    read_internal(hdu, I...)
read(hdu::ImageHDU, I::Int...) = read_internal(hdu, I...)[1]

read!(hdu::ImageHDU, array::StridedArray, I::Union{AbstractRange{Int}, Int, Colon}...) =
    read_internal!(hdu, array, I...)
read!(hdu::ImageHDU, array::StridedArray, I::Int...) = read_internal!(hdu, array, I...)[1]

"""
    write(f::FITS, data::StridedArray; header=nothing, name=nothing, ver=nothing)

Add a new image HDU to FITS file `f` with contents `data`. The
following array element types are supported: `UInt8`, `Int8`,
`UInt16`, `Int16`, `UInt32`, `Int32`, `Int64`, `Float32`,
`Float64`. If a `FITSHeader` object is passed as the `header` keyword
argument, the header will also be added to the new HDU.
"""
function write(f::FITS, data::StridedArray{T};
               header::Union{Nothing, FITSHeader}=nothing,
               name::Union{Nothing, String}=nothing,
               ver::Union{Nothing, Integer}=nothing) where {T<:Real}

    if !iscontiguous(data)
        throw(ArgumentError("data to be written out needs to be contiguously stored"))
    end

    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    if isa(header, FITSHeader)
        fits_write_header(f.fitsfile, header, true)
    end
    if isa(name, String)
        fits_update_key(f.fitsfile, "EXTNAME", name)
    end
    if isa(ver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", ver)
    end
    fits_write_pix(f.fitsfile, ones(Int, length(s)), length(data), data)
    nothing
end

"""
    write(hdu::ImageHDU, data::StridedArray)

Write data to an existing image HDU. 
The data to be written out must be stored contiguously in memory.
"""
function write(hdu::ImageHDU, data::StridedArray{<:Real})

    if !iscontiguous(data)
        throw(ArgumentError("data to be written out needs to be contiguously stored"))
    end

    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # Ensure sizes are equal
    hdu_size = fits_get_img_size(hdu.fitsfile)
    data_size = collect(size(data))

    if hdu_size != data_size
        error("size of HDU $(hdu_size) not equal to size of data $(data_size).")
    end

    fits_write_pix(hdu.fitsfile, ones(Int, length(size(data))), length(data), data)
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
    fits_assert_open(hdu.fitsfile)
    fits_assert_open(dest.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_copy_image_section(hdu.fitsfile, dest.fitsfile,
                            join([cfitsio_range_string(ri) for ri in r], ','))
end
