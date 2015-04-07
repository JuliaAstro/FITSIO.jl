import Base: checkbounds

# -----------------------------------------------------------------------------
# Types

# TODO : Cache metadata such as extname, extver, image size and data type?
#        This might allow faster access for size(ImageHDU) and ndim(ImageHDU)
abstract HDU
type ImageHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

type TableHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

type AsciiHDU <: HDU
    fitsfile::FITSFile
    ext::Int
end

# FITS is analagous to FITSFile, but holds a reference to all of its
# HDU objects. This is so that only a single HDU object is created for
# each extension in the file. It also allows a FITS object to tell
# previously created HDUs about events that happen to the file, such
# as deleting extensions. This could be done by, e.g., setting ext=-1 in
# the HDU object.
type FITS
    fitsfile::FITSFile
    filename::String
    mode::String
    hdus::Dict{Int, HDU}

    function FITS(filename::String, mode::String="r")
        f = (mode == "r"                     ? fits_open_file(filename, 0)   :
             mode == "r+" && isfile(filename)? fits_open_file(filename, 1)   :
             mode == "r+"                    ? fits_create_file(filename)    :
             mode == "w"                     ? fits_create_file("!"*filename):
             error("invalid open mode: $mode"))

        new(f, filename, mode, Dict{Int, HDU}())
    end
end

# FITSHeader stores the (key, value, comment) information for each card in
# a header. We could almost just use an OrderedDict for this, but we need
# to store comments so they are not lost when reading from disk,
# then writing back to disk with modifications.
type FITSHeader
    keys::Vector{ASCIIString}
    values::Vector{Any}
    comments::Vector{ASCIIString}
    map::Dict{ASCIIString, Int}

    function FITSHeader(keys::Vector{ASCIIString}, values::Vector{Any},
                        comments::Vector{ASCIIString})
        if ((length(keys) != length(values)) ||
            (length(keys) != length(comments)))
            error("keys, values, comments must be same length")
        end
        map = [keys[i]=>i for i=1:length(keys)]
        new(keys, values, comments, map)
    end
end

# -----------------------------------------------------------------------------
# FITS methods

function length(f::FITS)
    fits_assert_open(f.fitsfile)
    @compat Int(fits_get_num_hdus(f.fitsfile))
end

endof(f::FITS) = length(f)

function show(io::IO, f::FITS)
    fits_assert_open(f.fitsfile)

    print(io, "file: ", f.filename, "\n")
    print(io, "mode: ", f.mode, "\n")
    print(io, "extnum exttype         extname\n")

    for i = 1:length(f)
        hdutype = fits_movabs_hdu(f.fitsfile, i)
        extname = ""
        try
            extname = fits_read_keyword(f.fitsfile, "EXTNAME")[1]
        catch
            try
                extname = fits_read_keyword(f.fitsfile, "HDUNAME")[1]
            catch
            end
        end
        @printf io "%-6d %-15s %s\n" i hdutype extname
    end
end

# Returns HDU object based on extension number
function getindex(f::FITS, i::Integer)
    fits_assert_open(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    if i > length(f)
        error("index out of bounds")
    end
    hdutype = fits_movabs_hdu(f.fitsfile, i)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? AsciiHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end

# Returns HDU based on hduname, version
function getindex(f::FITS, name::String, ver::Int=0)
    fits_assert_open(f.fitsfile)
    fits_movnam_hdu(f.fitsfile, name, ver)
    i = fits_get_hdu_num(f.fitsfile)

    if haskey(f.hdus, i)
        return f.hdus[i]
    end

    hdutype = fits_get_hdu_type(f.fitsfile)
    f.hdus[i] = (hdutype == :image_hdu ? ImageHDU(f.fitsfile, i) :
                 hdutype == :binary_table ? TableHDU(f.fitsfile, i) :
                 hdutype == :ascii_table ? AsciiHDU(f.fitsfile, i) :
                 error("bad HDU type"))
    return f.hdus[i]
end


function close(f::FITS)
    fits_assert_open(f.fitsfile)
    fits_close_file(f.fitsfile)
    f.filename = ""
    f.mode = ""
    empty!(f.hdus)
    nothing
end

# -----------------------------------------------------------------------------
# Header methods

# returns one of: ASCIIString, Bool, Int, Float64, nothing
function parse_header_val(val::ASCIIString)
    len = length(val)
    if len < 1
        return nothing
    end
    c = val[1]
    if len == 1 && c == 'T'
        return true
    elseif len == 1 && c == 'F'
        return false
    elseif len >= 2 && c == '\'' && val[end] == '\''
        # This is a character string; according to FITS standard, trailing
        # spaces are insignificant, thus we remove them as well as the
        # surrounding quotes.
        return rstrip(val[2:end-1])
    else
        try
            return @compat parse(Int, val)
        catch
            try
                return @compat parse(Float64, val)
            catch
            end
        end
    end
    return val  # The value (probably) doesn't comply with the FITS
                # standard. Give up and return the unparsed string.
end

function readkey(fitsfile::FITSFile, key::Integer)
    fits_assert_open(fitsfile)
    keyout, value, comment = fits_read_keyn(fitsfile, key)
    keyout, parse_header_val(value), comment
end

function readkey(hdu::HDU, key::Integer)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    keyout, value, comment = fits_read_keyn(hdu.fitsfile, key)
    keyout, parse_header_val(value), comment
end

function readkey(fitsfile::FITSFile, key::ASCIIString)
    fits_assert_open(fitsfile)
    value, comment = fits_read_keyword(fitsfile, key)
    parse_header_val(value), comment
end

function readkey(hdu::HDU, key::ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    value, comment = fits_read_keyword(hdu.fitsfile, key)
    parse_header_val(value), comment
end

function readheader(fitsfile::FITSFile)
    fits_assert_open(fitsfile)

    # Below, we use a direct call to ffgkyn so that we can keep reusing the
    # same buffers.
    key = Array(Uint8, 81)
    value = Array(Uint8, 81)
    comment = Array(Uint8, 81)
    status = Cint[0]

    nkeys, morekeys = fits_get_hdrspace(fitsfile)

    # Initialize output arrays
    keys = Array(ASCIIString, nkeys)
    values = Array(Any, nkeys)
    comments = Array(ASCIIString, nkeys)
    for i=1:nkeys
        ccall((:ffgkyn,libcfitsio), Cint,
              (Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
              fitsfile.ptr, i, key, value, comment, status)
        keys[i] = bytestring(pointer(key))
        values[i] = parse_header_val(bytestring(pointer(value)))
        comments[i] = bytestring(pointer(comment))
    end
    fits_assert_ok(status[1])
    FITSHeader(keys, values, comments)
end

function readheader(hdu::HDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    readheader(hdu.fitsfile)
end

length(hdr::FITSHeader) = length(hdr.keys)
haskey(hdr::FITSHeader, key::ASCIIString) = in(key, hdr.keys)
keys(hdr::FITSHeader) = hdr.keys
values(hdr::FITSHeader) = hdr.values
getindex(hdr::FITSHeader, key::ASCIIString) = hdr.values[hdr.map[key]]
getindex(hdr::FITSHeader, i::Integer) = hdr.values[i]

function setindex!(hdr::FITSHeader, value::Any, key::ASCIIString)
    if in(key, hdr.keys)
        hdr.values[hdr.map[key]] = value
    else
        push!(hdr.keys, key)
        push!(hdr.values, value)
        push!(hdr.comments, "")
        hdr.map[key] = length(hdr.keys)
    end
end

function setindex!(hdr::FITSHeader, value::Any, i::Integer)
    hdr.values[i] = value
end

# Comments
getcomment(hdr::FITSHeader, key::ASCIIString) = hdr.comments[hdr.map[key]]
getcomment(hdr::FITSHeader, i::Integer) = hdr.comments[i]
function setcomment!(hdr::FITSHeader, key::ASCIIString, comment::ASCIIString)
    hdr.comments[hdr.map[key]] = comment
end
function setcomment!(hdr::FITSHeader, i::Integer, comment::ASCIIString)
    hdr.comments[i] = comment
end

# Display the header
hdrval_to_str(val::Bool) = val ? "T" : "F"
hdrval_to_str(val::ASCIIString) = @sprintf "'%s'" val
hdrval_to_str(val::Union(FloatingPoint, Integer)) = string(val)
function show(io::IO, hdr::FITSHeader)
    for i=1:length(hdr)
        cl = length(hdr.comments[i])
        if hdr.keys[i] == "COMMENT" || hdr.keys[i] == "HISTORY"
            if cl > 71
                @printf io "%s %s\n" hdr.keys[i] hdr.comments[i][1:71]
            else
                @printf io "%s %s\n" hdr.keys[i] hdr.comments[i]
            end
        else
            @printf io "%-8s" hdr.keys[i]
            if hdr.values[i] === nothing
                print(io, "                      ")
                rc = 50  # remaining characters on line
            else
                val = hdrval_to_str(hdr.values[i])
                @printf io "= %20s" val
                rc = length(val) <= 20 ? 50 : 70 - length(val)
            end

            if cl > 0
                if cl > rc - 3
                    @printf io " / %s" hdr.comments[i][1:rc-3]
                else
                    @printf io " / %s" hdr.comments[i]
                end
            end
            print(io, "\n")
        end
    end
end

const RESERVED_KEYS = ["SIMPLE","EXTEND","XTENSION","BITPIX","PCOUNT","GCOUNT",
                       "THEAP","EXTNAME","BUNIT","BSCALE","BZERO","BLANK",
                       "ZQUANTIZ","ZDITHER0","ZIMAGE","ZCMPTYPE","ZSIMPLE",
                       "ZTENSION","ZPCOUNT","ZGCOUNT","ZBITPIX","ZEXTEND",
                       "CHECKSUM","DATASUM"]

# This is more complex than you would think because some reserved keys
# are only reserved when other keys are present. Also, in general a key
# may appear more than once in a header.
function reserved_key_indicies(hdr::FITSHeader)
    nhdr = length(hdr)
    indicies = Int[]
    for i=1:nhdr
        if in(hdr.keys[i], RESERVED_KEYS)
            push!(indicies, i)
        end
    end

    # Note that this removes anything matching NAXIS\d regardless of # of axes.
    if in("NAXIS", hdr.keys)
        for i=1:nhdr
            if ismatch(r"^NAXIS\d*$", hdr.keys[i])
                push!(indicies, i)
            end
        end
    end

    if in("ZNAXIS", hdr.keys)
        for i=1:nhdr
            if (ismatch(r"^ZNAXIS\d*$", hdr.keys[i]) ||
                ismatch(r"^ZTILE\d*$", hdr.keys[i]) ||
                ismatch(r"^ZNAME\d*$", hdr.keys[i]) ||
                ismatch(r"^ZVAL\d*$", hdr.keys[i]))
                push!(indicies, i)
            end
        end
    end

    if in("TFIELDS", hdr.keys)
        for i=1:nhdr
            for re in [r"^TFORM\d*$", r"^TTYPE\d*$", r"^TDIM\d*$",
                       r"^TUNIT\d*$", r"^TSCAL\d*$", r"^TZERO\d*$",
                       r"^TNULL\d*$", r"^TDISP\d*$", r"^TDMIN\d*$",
                       r"^TDMAX\d*$", r"^TDESC\d*$", r"^TROTA\d*$",
                       r"^TRPIX\d*$", r"^TRVAL\d*$", r"^TDELT\d*$",
                       r"^TCUNI\d*$", r"^TFIELDS$"]
                if ismatch(re, hdr.keys[i])
                    push!(indicies, i)
                end
            end
        end
    end

    return indicies
end


# Header writing: "low-level" (works directly on FITSFile)
# if `clean` is true, skip writing reserved header keywords
function write_header(f::FITSFile, hdr::FITSHeader, clean::Bool=true)
    indicies = clean? reserved_key_indicies(hdr): Int[]
    for i=1:length(hdr)
        if clean && in(i, indicies)
            continue
        end
        if hdr.keys[i] == "COMMENT"
            fits_write_comment(f, hdr.comments[i])
        elseif hdr.keys[i] == "HISTORY"
            fits_write_history(f, hdr.comments[i])
        elseif hdr.comments[i] == ""
            fits_update_key(f, hdr.keys[i], hdr.values[i])
        else
            fits_update_key(f, hdr.keys[i], hdr.values[i], hdr.comments[i])
        end
    end
end

# -----------------------------------------------------------------------------
# ImageHDU methods

# Display the image datatype and dimensions
function show(io::IO, hdu::ImageHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    sz = fits_get_img_size(hdu.fitsfile)
    @printf io "file: %s\nextension: %d\ntype: IMAGE\nimage info:\n  bitpix: %d\n  size: %s" fits_file_name(hdu.fitsfile) hdu.ext bitpix tuple(sz...)
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
    data = Array(bitpix_to_type[bitpix], sz...)
    fits_read_pix(hdu.fitsfile, data)
    data
end

# `trailingsize` is used to compute last index in checkbounds.
# sz: array size (tuple), n: starting index.
# Same as Base.trailingsize but takes size tuple rather than array.
function trailingsize(sz, n)
    s = 1
    for i=n:length(sz)
        s *= sz[i]
    end
    return s
end

# `checkbounds` method, used for ImageHDU. (ImageHDU is not a subtype
# of AbstractArray, so we can't just use methods in Base.) Note that
# this takes a size tuple whereas Methods in Base take array A and use
# size(A,n). For ImageHDU, size(hdu,n) gets size(hdu) as in
# intermediary, so it is better to just get that tuple once.
function checkbounds(sz::NTuple, I::Union(Int, Range{Int})...)
    n = length(I)
    if n > 0
        for dim = 1:(n-1)
            checkbounds(sz[dim], I[dim])
        end
        checkbounds(trailingsize(sz,n), I[n])
    end
end

# Read a subset of an ImageHDU
function _getindex(hdu::ImageHDU, I::Union(Range{Int},Int)...)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # check number of indicies and bounds.
    # Note that N_indicies != ndim is not supported, differing from Array.
    # It could be supported in the future with care taken in constructing
    # first, last, step arrays passed to cfitsio.
    sz = tuple(fits_get_img_size(hdu.fitsfile)...)
    if length(I) != length(sz)
        throw(DimensionMismatch("number of indicies must match dimensions"))
    end
    checkbounds(sz, I...)

    # construct first, last and step vectors
    firsts = Clong[first(i) for i in I]
    lasts = Clong[last(i) for i in I]
    steps = Clong[isa(i, Int)? 1: step(i) for i in I]

    # construct output array
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    data = Array(bitpix_to_type[bitpix], Base.index_shape(I...))

    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

# general method and version that returns a single value rather than 0-d array
getindex(hdu::ImageHDU, I::Union(Range{Int},Int)...) = _getindex(hdu, I...)
getindex(hdu::ImageHDU, I::Int...) = _getindex(hdu, I...)[1]

# Add a new ImageHDU to a FITS object
# The following Julia data types are supported for writing images by cfitsio:
# Uint8, Int8, Uint16, Int16, Uint32, Int32, Int64, Float32, Float64
function write{T}(f::FITS, data::Array{T};
                  header::Union(Nothing, FITSHeader)=nothing)
    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    if isa(header, FITSHeader)
        write_header(f.fitsfile, header, true)
    end
    fits_write_pix(f.fitsfile, ones(Int, length(s)), length(data), data)
    nothing
end

# Copy a rectangular section of an image and write it to a new FITS
# primary image or image extension. The new image HDU is appended to
# the end of the output file; all the keywords in the input image will
# be copied to the output image. The common WCS keywords will be
# updated if necessary to correspond to the coordinates of the section.

# TODO: Change Range types once v0.2 is no longer supported.

range2fits_str(r::Range1) = @sprintf "%d:%d" first(r) last(r)
range2fits_str(r::Range) = @sprintf "%d:%d:%d" first(r) last(r) step(r)
fits_copy_image_section(fin::FITSFile, fout::FITSFile, r...) =
    fits_copy_image_section(fin, fout, join([range2str(ri) for ri in r], ','))
function copy_section(hdu::ImageHDU, destination::FITS, r::Range...)
    fits_assert_open(hdu.fitsfile)
    fits_assert_open(destination.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_copy_image_section(hdu.fitsfile, destination.fitsfile,
                            join([range2fits_str(ri) for ri in r], ','))
end
