using .Libcfitsio

# There are a few direct `ccall`s to libcfitsio in this module. For this, we
# need a few non-exported things from Libcfitsio: the shared library handle,
# and a helper function for raising errors. TYPE_FROM_BITPIX is awkwardly
# defined in Libcfitsio, even though it is not used there.
import .Libcfitsio: libcfitsio,
                    fits_assert_ok,
                    TYPE_FROM_BITPIX

function libcfitsio_version()
    # fits_get_version returns a float. e.g., 3.341f0. We parse that
    # into a proper version number. E.g., 3.341 -> v"3.34.1"
    v = convert(Int, round(1000 * fits_get_version()))
    x = div(v, 1000)
    y = div(rem(v, 1000), 10)
    z = rem(v, 10)
    VersionNumber(x, y, z)
end

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

type ASCIITableHDU <: HDU
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

    function FITSHeader(keys::Vector{ASCIIString}, values::Vector,
                        comments::Vector{ASCIIString})
        if ((length(keys) != length(values)) ||
            (length(keys) != length(comments)))
            error("keys, values, comments must be same length")
        end
        map = [keys[i]=>i for i=1:length(keys)]
        new(keys, convert(Vector{Any}, values), comments, map)
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

    # Get name and type of all HDUs.
    nhdu = length(f)
    names = Array(ASCIIString, nhdu)
    vers = Array(ASCIIString, nhdu)
    types = Array(ASCIIString, nhdu)
    for i = 1:nhdu
        t = fits_movabs_hdu(f.fitsfile, i)
        types[i] = (t == :image_hdu ? "Image" :
                    t == :binary_table ? "Table" :
                    t == :ascii_table ? "ASCIITable" :
                    error("unknown HDU type"))
        nname = _try_read_key(f.fitsfile, ASCIIString, ("EXTNAME", "HDUNAME"))
        names[i] = get(nname, "")
        nver = _try_read_key(f.fitsfile, Int, ("EXTVER", "HDUVER"))
        vers[i] = isnull(nver) ? "" : string(get(nver))
    end

    namelen = max(maximum(length, names), 4) + 2
    verlen = maximum(length, vers)
    if verlen > 0
        verlen = max(verlen, 3) + 2
        namehead = string(rpad("Name", namelen), rpad("Ver", verlen))
    else
        namehead = rpad("Name", namelen)
    end

    print(io, """File: $(f.filename)
    Mode: \"$(f.mode)\"
    HDUs: Num   $(namehead)Type
    """)
    for i in 1:nhdu
        @printf io "      %-5d %s%s%s\n" i rpad(names[i], namelen) rpad(vers[i], verlen) types[i]
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
                 hdutype == :ascii_table ? ASCIITableHDU(f.fitsfile, i) :
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
                 hdutype == :ascii_table ? ASCIITableHDU(f.fitsfile, i) :
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

# modified version of Base.rstrip that returns " " rather than "" if the
# string is all spaces. This is necessary because trailing whitespace is
# not significant in FITS header keywords, but *leading* whitespace is.
# See CFITSIO manual section 4.5 for details.
function _rstrip_fits(s::ASCIIString)

    # return input string if length is 0 or 1 or last char is not space.
    if length(s) == 0 || length(s) == 1 || s[end] != ' '
        return s
    end

    i = endof(s) - 1
    while i > 1
        if s[i] != ' '
            return s[1:i]
        end
        i -= 1
    end
    return s[1:1]
end

# parse FITS header value as an ASCIIString.
function parse_header_val(::Type{ASCIIString}, val::ASCIIString)
    if length(val) < 2 || val[1] != '\'' || val[end] != '\''
        throw(ArgumentError("ASCIIString header values must begin and end with ' character"))
    end
    return _rstrip_fits(val[2:end-1])
end

parse_header_val(::Type{Int}, val::ASCIIString) = parse(Int, val)

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
        return _rstrip_fits(val[2:end-1])
    else
        # TODO: change this block to `tryparse` once support for
        # tryparse(Int, x) is in Compat or v0.3 is no longer supported.
        try
            return parse(Int, val)
        catch
        end

        fval = tryparse(Float64, val)
        isnull(fval) || return get(fval)
    end

    # The value (probably) doesn't comply with the FITS standard.
    # Give up and return the unparsed string.
    return val

end

function read_key(hdu::HDU, key::Integer)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    keyout, value, comment = fits_read_keyn(hdu.fitsfile, key)
    keyout, parse_header_val(value), comment
end

function read_key(hdu::HDU, key::ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    value, comment = fits_read_keyword(hdu.fitsfile, key)
    parse_header_val(value), comment
end

# Try to read the raw keys in order given; returns Nullable{ASCIIString}.
# (null if key doesn't exist.)
function _try_read_key{T}(f::FITSFile, ::Type{T}, names)
    status = Cint[0]
    value = Array(Uint8, 71)
    for name in names
        ccall((:ffgkey, libcfitsio), Cint,
              (Ptr{Void},Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
              f.ptr, bytestring(name), value, C_NULL, status)

        # If the key is found, return it. If there was some other error
        # besides key not found, throw an error.
        if status[1] == 0
            return Nullable(parse_header_val(T, bytestring(pointer(value))))
        elseif status[1] != 202
            error(fits_get_errstatus(status[1]))
        end
    end
    return Nullable{ASCIIString}()
end

function read_header(hdu::HDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # Below, we use a direct call to ffgkyn so that we can keep reusing the
    # same buffers.
    key = Array(Uint8, 81)
    value = Array(Uint8, 81)
    comment = Array(Uint8, 81)
    status = Cint[0]

    nkeys, morekeys = fits_get_hdrspace(hdu.fitsfile)

    # Initialize output arrays
    keys = Array(ASCIIString, nkeys)
    values = Array(Any, nkeys)
    comments = Array(ASCIIString, nkeys)
    for i=1:nkeys
        ccall((:ffgkyn,libcfitsio), Cint,
              (Ptr{Void},Cint,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Cint}),
              hdu.fitsfile.ptr, i, key, value, comment, status)
        keys[i] = bytestring(pointer(key))
        values[i] = parse_header_val(bytestring(pointer(value)))
        comments[i] = bytestring(pointer(comment))
    end
    fits_assert_ok(status[1])
    FITSHeader(keys, values, comments)
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
get_comment(hdr::FITSHeader, key::ASCIIString) = hdr.comments[hdr.map[key]]
get_comment(hdr::FITSHeader, i::Integer) = hdr.comments[i]
function set_comment!(hdr::FITSHeader, key::ASCIIString, comment::ASCIIString)
    hdr.comments[hdr.map[key]] = comment
end
function set_comment!(hdr::FITSHeader, i::Integer, comment::ASCIIString)
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

# helper function for show(::HDU) - gets string with HDU number, name, version.
# assumes file is open or and on correct current hdu.
function _get_hdu_info_string(hdu::HDU)
    hduname = _try_read_key(hdu.fitsfile, ASCIIString, ("EXTNAME", "HDUNAME"))
    hduver = _try_read_key(hdu.fitsfile, Int, ("EXTVER", "HDUVER"))
    if !isnull(hduname) && !isnull(hduver)
        return "$(hdu.ext) (name=$(repr(get(hduname))), ver=$(get(hduver)))"
    elseif !isnull(hduname)
        return "$(hdu.ext) (name=$(repr(get(hduname))))"
    end
    return string(hdu.ext)
end

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
        datainfo = @sprintf "%s (physical: %s)" TYPE_FROM_BITPIX[equivbitpix] TYPE_FROM_BITPIX[bitpix]
    end
    
    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(_get_hdu_info_string(hdu))
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
function _read(hdu::ImageHDU, I::Union(Range{Int},Int, Colon)...)
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
read(hdu::ImageHDU, I::Union(Range{Int}, Int, Colon)...) = _read(hdu, I...)
read(hdu::ImageHDU, I::Int...) = _read(hdu, I...)[1]

# Add a new ImageHDU to a FITS object
# The following Julia data types are supported for writing images by cfitsio:
# Uint8, Int8, Uint16, Int16, Uint32, Int32, Int64, Float32, Float64
function write{T}(f::FITS, data::Array{T};
                  header::Union(Nothing, FITSHeader)=nothing,
                  hduname::Union(Nothing, ASCIIString)=nothing,
                  hduver::Union(Nothing, Integer)=nothing)
    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    if isa(header, FITSHeader)
        write_header(f.fitsfile, header, true)
    end
    if isa(hduname, ASCIIString)
        fits_update_key(f.fitsfile, "EXTNAME", hduname)
    end
    if isa(hduver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", hduver)
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
cfitsio_string(r::UnitRange) = @sprintf "%d:%d" first(r) last(r)
cfitsio_string(r::StepRange) = @sprintf "%d:%d:%d" first(r) last(r) step(r)

function copy_section(hdu::ImageHDU, dest::FITS, r::Range{Int}...)
    fits_assert_open(hdu.fitsfile)
    fits_assert_open(dest.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_copy_image_section(hdu.fitsfile, dest.fitsfile,
                            join([cfitsio_string(ri) for ri in r], ','))
end


# -----------------------------------------------------------------------------
# TableHDU methods

# Table type conversions:
#
# fits_tform_letter(::DataType) -> Char
#     Given Julia array eltype, what FITS table letter code should be used
#     when defining a table column? For example, to store an array of UInt32,
#     use the TFORM letter 'V'.
# CFITSIO_COLTYPE[::Int] -> DataType
#     Given return code from fits_get_eqcoltype(), what type of Julia array
#     should be constructed? For example, for 'V' columns,
#     fits_get_eqcoltype() returns 40. This function maps that code back to
#     UInt32. This also illustrates why we can't simply use the normal CFITSIO
#     datatype mapping: 40 would map to Culong, which is a 64-bit unsigned
#     integer on 64-bit UNIX platforms.
const CFITSIO_COLTYPE = Dict{Int, DataType}()
for (T, tform, code) in ((UInt8,       'B',  11),
                         (Int8,        'S',  12),
                         (Bool,        'L',  14),
                         (ASCIIString, 'A',  16),
                         (UInt16,      'U',  20),
                         (Int16,       'I',  21),
                         (UInt32,      'V',  40),
                         (Int32,       'J',  41),
                         (Int64,       'K',  81),
                         (Float32,     'E',  42),
                         (Float64,     'D',  82),
                         (Complex64,   'C',  83),
                         (Complex128,  'M', 163))
    @eval fits_tform_char(::Type{$T}) = $tform
    CFITSIO_COLTYPE[code] = T
end
typealias FITSTableScalar Union(UInt8, Int8, Bool, UInt16, Int16, Uint32,
                                Int32, Int64, Float32, Float64, Complex64,
                                Complex128)

# Helper function for reading information about a (binary) table column
# Returns: (eltype, rowsize, isvariable)
function _get_col_info(f::FITSFile, colnum::Integer)
    eqtypecode, repeat, width = fits_get_eqcoltype(f, colnum)
    isvariable = eqtypecode < 0
    eqtypecode = abs(eqtypecode)

    (eqtypecode == 1) && error("BitArray ('X') columns not yet supported")

    T = CFITSIO_COLTYPE[eqtypecode]

    if isvariable
        if T !== ASCIIString
            T = Vector{T}
        end
        rowsize = Int[]
    else
        if T === ASCIIString
            # for strings, cfitsio only considers it to be a vector column if
            # width != repeat, even if tdim is multi-valued.
            if repeat == width
                rowsize = Int[]
            else
                tdim = fits_read_tdim(f, colnum)
                # if tdim isn't multi-valued, ignore it (we know it *is* a
                # vector column). If it is multi-valued, prefer it to repeat
                # width.
                if length(tdim) == 1
                    rowsize = [div(repeat, width)]
                else
                    rowsize = tdim[2:end]
                end
            end
        else
            if repeat == 1
                rowsize = Int[]
            else
                rowsize = fits_read_tdim(f, colnum)
            end
        end
    end

    return T, rowsize, isvariable
end

# Helper function for getting fits tdim shape for given array
fits_tdim(A::Array) = (ndims(A) == 1)? [1]: [size(A, i) for i=1:ndims(A)-1]
function fits_tdim(A::Array{ASCIIString})
    n = ndims(A)
    tdim = Array(Int, n)
    tdim[1] = maximum(length, A)
    for i=2:n
        tdim[n] = size(A, n-1)
    end
    tdim
end

# Helper function for getting fits tform string for given table type
# and data array.
fits_tform{T}(::Type{TableHDU}, A::Array{T}) = "$(prod(fits_tdim(A)))$(fits_tform_char(T))"

# For string arrays with 2+ dimensions, write tform as rAw. Otherwise,
# cfitsio doesn't recognize that multiple strings should be written to
# a single row, even if TDIM is set to 2+ dimensions.
fits_tform(::Type{TableHDU}, A::Vector{ASCIIString}) = "$(maximum(length, A))A"
fits_tform(::Type{TableHDU}, A::Array{ASCIIString}) = "$(prod(fits_tdim(A)))A$(maximum(length, A))"

# variable length columns
fits_tform_v{T<:FITSTableScalar}(::Type{TableHDU}, A::Vector{Vector{T}}) = "1P$(fits_tform_char(T))($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{ASCIIString}) = "1PA($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{Vector}) = error("column data must be a leaf type: e.g., Vector{Vector{Int}}, not Vector{Vector{T}}.")
fits_tform_v(::Type{TableHDU}, ::Any) = error("variable length columns only supported for arrays of arrays and arrays of ASCIIString")
fits_tform_v(::Type{ASCIITableHDU}, ::Any) = error("variable length columns not supported in ASCII tables")


fits_tform(::Type{ASCIITableHDU}, ::Vector{Int16}) = "I7"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Int32}) = "I12"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float32}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float64}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, A::Vector{ASCIIString}) = "A$(maximum(length, A))"
fits_tform(::Type{ASCIITableHDU}, A::Vector) = error("unsupported type: $(eltype(A))")
fits_tform(::Type{ASCIITableHDU}, A::Array) = error("only 1-d arrays supported: dimensions are $(size(A))")

# for passing to fits_create_tbl.
table_type_code(::Type{ASCIITableHDU}) = convert(Cint, 1)
table_type_code(::Type{TableHDU}) = convert(Cint, 2)

function show(io::IO, hdu::TableHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    # allocate return arrays for column names & types
    colnames_in = [Array(Uint8, 70) for i=1:ncols]
    coltypes_in = [Array(Uint8, 70) for i=1:ncols]
    nrows = Array(Int64, 1)
    status = Cint[0]

    # fits_read_btblhdrll (Can pass NULL for return fields not needed.)
    ccall(("ffghbnll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Cint}, Ptr{Ptr{Uint8}},  # nrows, tfields, ttype
           Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Ptr{Uint8},  # tform,tunit,extname
           Ptr{Clong}, Ptr{Cint}),  # pcount, status
          hdu.fitsfile.ptr, ncols, nrows, C_NULL, colnames_in, coltypes_in,
          C_NULL, C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    colnames = [bytestring(pointer(item)) for item in colnames_in]
    coltypes = [bytestring(pointer(item)) for item in coltypes_in]

    maxlen = maximum(length, colnames)

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(_get_hdu_info_string(hdu))
    Type: Table
    Rows: $(nrows[1])
    Columns: """)
    print(io, rpad("Name", maxlen), "  Type\n")
    for i in 1:ncols
        T, rowsize, isvariable = _get_col_info(hdu.fitsfile, i)

        print(io, "         $(rpad(colnames[i], maxlen))  $T")
        if length(rowsize) > 0
            print(io, "  $(tuple(rowsize...))")
        end
        if isvariable
            print(io, "  [Var]")
        end
        print(io, "\n")
    end
end

function show(io::IO, hdu::ASCIITableHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    # allocate return arrays for column names & types
    colnames_in = [Array(Uint8, 70) for i=1:ncols]
    coltypes_in = [Array(Uint8, 70) for i=1:ncols]
    nrows = Array(Int64, 1)
    status = Cint[0]

    # fits_read_atblhdrll (Can pass NULL for return fields not needed)
    ccall(("ffghtbll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Int64}, Ptr{Cint},  # rowlen, nrows, tfields
           Ptr{Ptr{Uint8}}, Ptr{Clong}, Ptr{Ptr{Uint8}},  # ttype, tbcol, tform
           Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Cint}),  # tunit, extname, status
          hdu.fitsfile.ptr, ncols,
          C_NULL, nrows, C_NULL,
          colnames_in, C_NULL, coltypes_in,
          C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    colnames = [bytestring(pointer(item)) for item in colnames_in]
    coltypes = [bytestring(pointer(item)) for item in coltypes_in]

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(_get_hdu_info_string(hdu))
    Type: ASCII TABLE
    Rows: $(nrows[1])
    Columns:
    """)
    for i in 1:ncols
        @printf io "    %s  (%s)\n" colnames[i] coltypes[i]
    end
end

# Write a variable length array column of numbers
# (separate implementation from normal fits_write_col function because
#  we must make separate calls to `fits_write_col` for each row.)
function _write_var_col{T}(f::FITSFile, colnum::Integer,
                           data::Vector{Vector{T}})
    for i=1:length(data)
        fits_write_col(f, colnum, i, 1, data[i])
    end
end

function _write_var_col(f::FITSFile, colnum::Integer,
                        data::Vector{ASCIIString})
    status = Cint[0]
    buffer = Array(Ptr{UInt8}, 1)  # holds the address of the current row
    for i=1:length(data)
        buffer[1] = pointer(data[i])

        # Note that when writing to a variable ASCII column, the
        # ‘firstelem’ and ‘nelements’ parameter values in the
        # fits_write_col routine are ignored and the number of
        # characters to write is simply determined by the length of
        # the input null-terminated character string.
        ccall((:ffpcls, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Int64, Int64, Ptr{Ptr{Uint8}},
               Ptr{Cint}),
              f.ptr, colnum, i, 1, length(data[i]), buffer, status)
        fits_assert_ok(status[1])
    end
end

# Add a new TableHDU to a FITS object
function write_impl(f::FITS, colnames::Vector{ASCIIString}, coldata::Vector,
                    hdutype, hduname, hduver, header, units, varcols)
    fits_assert_open(f.fitsfile)

    # move to last HDU; table will be added after the CHDU
    nhdus = @compat(Int(fits_get_num_hdus(f.fitsfile)))
    (nhdus > 1) && fits_movabs_hdu(f.fitsfile, nhdus)

    ncols = length(colnames)
    ttype = [pointer(name) for name in colnames]

    # determine which columns are requested to be variable-length
    isvarcol = zeros(Bool, ncols)
    if !isa(varcols, Nothing)
        for i=1:ncols
            isvarcol[i] = (i in varcols) || (colnames[i] in varcols)
        end
    end

    # create an array of tform strings (which we will create pointers to)
    tform_str = Array(ASCIIString, ncols)
    for i in 1:ncols
        if isvarcol[i]
            tform_str[i] = fits_tform_v(hdutype, coldata[i])
        else
            tform_str[i] = fits_tform(hdutype, coldata[i])
        end
    end
    tform = [pointer(s) for s in tform_str]

    # get units
    if isa(units, Nothing)
        tunit = C_NULL
    else
        tunit = Ptr{Uint8}[(haskey(units, name)? pointer(units[name]): C_NULL)
                           for name in colnames]
    end

    # extension name
    hduname_ptr = (isa(hduname, Nothing) ? convert(Ptr{Uint8}, C_NULL) :
                   pointer(hduname))

    status = Cint[0]
    ccall(("ffcrtb", libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Cint, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}},
           Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Cint}),
          f.fitsfile.ptr, table_type_code(hdutype), 0, ncols,  # 0 = nrows
          ttype, tform, tunit, hduname_ptr, status)
    fits_assert_ok(status[1])

    # For binary tables, write tdim info
    if hdutype === TableHDU
        for (i, a) in enumerate(coldata)
            isvarcol[i] || fits_write_tdim(f.fitsfile, i, fits_tdim(a))
        end
    end

    if isa(header, FITSHeader)
        write_header(f.fitsfile, header, true)
    end
    if isa(hduver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", hduver)
    end

    for (i, a) in enumerate(coldata)
        if isvarcol[i]
            _write_var_col(f.fitsfile, i, a)
        else
            fits_write_col(f.fitsfile, i, 1, 1, a)
        end
    end
    nothing
end

function write(f::FITS, colnames::Vector{ASCIIString}, coldata::Vector;
               units=nothing, header=nothing, hdutype=TableHDU,
               hduname=nothing, hduver=nothing, varcols=nothing)
    if length(colnames) != length(coldata)
        error("length of colnames and length of coldata must match")
    end
    write_impl(f, colnames, coldata, hdutype, hduname, hduver, header,
               units, varcols)
end

function write(f::FITS, data::Dict{ASCIIString};
               units=nothing, header=nothing, hdutype=TableHDU,
               hduname=nothing, hduver=nothing, varcols=nothing)
    colnames = collect(keys(data))
    coldata = collect(values(data))
    write_impl(f, colnames, coldata, hdutype, hduname, hduver, header,
               units, varcols)
end

# Read a variable length array column of numbers
# (separate implementation from normal fits_read_col function because
# the length of each vector must be determined for each row.
function _read_var_col{T}(f::FITSFile, colnum::Integer, data::Vector{Vector{T}})
    nrows = length(data)
    for i=1:nrows
        repeat, offset = fits_read_descript(f, colnum, i)
        data[i] = Array(T, repeat)
        fits_read_col(f, colnum, i, 1, data[i])
    end
end

# Read a variable length array column of strings
# (Must be separate implementation from normal fits_read_col function because
# the length of each string must be determined for each row.)
function _read_var_col(f::FITSFile, colnum::Integer, data::Vector{ASCIIString})
    status = Cint[0]
    bufptr = Array(Ptr{UInt8}, 1)  # holds a pointer to the current row buffer
    for i=1:length(data)
        repeat, offset = fits_read_descript(f, colnum, i)
        buffer = Array(UInt8, repeat)
        bufptr[1] = pointer(buffer)
        ccall((:ffgcvs, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Int64, Int64,
               Ptr{Uint8}, Ptr{Ptr{Uint8}}, Ptr{Cint}, Ptr{Cint}),
              f.ptr, colnum, i, 1, repeat, " ", bufptr, C_NULL, status)
        fits_assert_ok(status[1])

        # Create string out of the buffer, terminating at null characters
        zeropos = search(buffer, 0x00)
        data[i] = (zeropos >= 1) ? ASCIIString(buffer[1:(zeropos-1)]) :
                                   ASCIIString(buffer)
    end
end

# Read a table column
function read(hdu::ASCIITableHDU, colname::ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname)

    # `eqcoltype`: do SCALE/ZERO conversion automatically
    typecode, repcnt, width = fits_get_eqcoltype(hdu.fitsfile, colnum)

    T = CFITSIO_COLTYPE[typecode]
    result = Array(T, nrows)
    fits_read_col(hdu.fitsfile, colnum, 1, 1, result)

    return result
end

function read(hdu::TableHDU, colname::ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname)

    T, rowsize, isvariable = _get_col_info(hdu.fitsfile, colnum)

    result = Array(T, rowsize..., nrows)

    if isvariable
        _read_var_col(hdu.fitsfile, colnum, result)
    else
        fits_read_col(hdu.fitsfile, colnum, 1, 1, result)
    end

    return result
end
