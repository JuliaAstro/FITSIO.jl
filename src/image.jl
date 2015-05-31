# ImageHDU methods

# Display the image datatype and dimensions
function show(io::IO, hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_type(hdu.fitsfile)
    equivbitpix = fits_get_img_equivtype(hdu.fitsfile)
    sz = fits_get_img_size(hdu.fitsfile)

    if bitpix == equivbitpix
        datainfo = string(TYPE_FROM_BITPIX[equivbitpix])
    else
        datainfo = @sprintf("%s (physical: %s)",
                            TYPE_FROM_BITPIX[equivbitpix],
                            TYPE_FROM_BITPIX[bitpix])
    end
    
    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(hdu.ext)$(fits_get_ext_info_string(hdu.fitsfile))
    Type: Image
    Datatype: $datainfo
    Datasize: $(tuple(sz...))
    """)
end

# Get image dimensions
function ndims(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_get_img_dim(hdu.fitsfile)
end

# Get image size
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

# `endof` is needed so that hdu[:] can throw DimensionMismatch
# when ndim != 1, rather than no method.
length(hdu::ImageHDU) = prod(size(hdu))
endof(hdu::ImageHDU) = length(hdu::ImageHDU)

# Read a full image from an HDU
function read(hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    sz = fits_get_img_size(hdu.fitsfile)
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    data = Array(TYPE_FROM_BITPIX[bitpix], sz...)
    fits_read_pix(hdu.fitsfile, data)
    data
end

# `trailingsize` is used to compute last index for checkbounds.
# sz: array size (tuple), n: starting index.
# Same as Base.trailingsize but takes size tuple rather than array.
function trailingsize(sz, n)
    s = 1
    for i=n:length(sz)
        s *= sz[i]
    end
    return s
end

# Read a subset of an ImageHDU
function read_internal(hdu::ImageHDU, I::Union(Range{Int},Int, Colon)...)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # check number of indicies and bounds.
    # Note that N_indicies != ndim is not supported, differing from Array.
    # It could be supported in the future with care taken in constructing
    # first, last, step arrays passed to cfitsio.
    sz = tuple(fits_get_img_size(hdu.fitsfile)...)
    n = length(I)
    if n != length(sz)
        throw(DimensionMismatch("number of indicies must match dimensions"))
    end

    # colon-expanded version of I
    J = ntuple(n, i-> isa(I[i], Colon)? (1:sz[i]) : I[i])

    if n > 0
        for i = 1:(n-1)
            Base.checkbounds(sz[i], J[i])
        end
        Base.checkbounds(trailingsize(sz,n), J[n])
    end

    # construct first, last and step vectors
    firsts = Clong[first(idx) for idx in J]
    lasts = Clong[last(idx) for idx in J]
    steps = Clong[isa(idx, Range)? step(idx): 1 for idx in J]

    # construct output array
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    data = Array(TYPE_FROM_BITPIX[bitpix], Base.index_shape(J...))

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

# general method and version that returns a single value rather than 0-d array
read(hdu::ImageHDU, I::Union(Range{Int}, Int, Colon)...) =
    read_internal(hdu, I...)
read(hdu::ImageHDU, I::Int...) = read_internal(hdu, I...)[1]

# Add a new ImageHDU to a FITS object
# The following Julia data types are supported for writing images by cfitsio:
# Uint8, Int8, Uint16, Int16, Uint32, Int32, Int64, Float32, Float64
function write{T}(f::FITS, data::Array{T};
                  header::Union(Nothing, FITSHeader)=nothing,
                  name::Union(Nothing, ASCIIString)=nothing,
                  ver::Union(Nothing, Integer)=nothing)
    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    if isa(header, FITSHeader)
        fits_write_header(f.fitsfile, header, true)
    end
    if isa(name, ASCIIString)
        fits_update_key(f.fitsfile, "EXTNAME", name)
    end
    if isa(ver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", ver)
    end
    fits_write_pix(f.fitsfile, ones(Int, length(s)), length(data), data)
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

function copy_section(hdu::ImageHDU, dest::FITS, r::Range{Int}...)
    fits_assert_open(hdu.fitsfile)
    fits_assert_open(dest.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_copy_image_section(hdu.fitsfile, dest.fitsfile,
                            join([cfitsio_range_string(ri) for ri in r], ','))
end
