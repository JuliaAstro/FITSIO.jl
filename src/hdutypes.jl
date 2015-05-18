using .Libcfitsio

# There are a few direct `ccall`s to libcfitsio in this module. For this, we
# need a few non-exported things from Libcfitsio: the shared library handle,
# and a helper function for raising errors. TYPE_FROM_BITPIX is awkwardly
# defined in Libcfitsio, even though it is not used there.
import .Libcfitsio: libcfitsio,
                    fits_assert_ok,
                    TYPE_FROM_BITPIX

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
    CFITSIO_COLTYPE[-code] = T  # variable length arrays
end

## Helper functions for writing a table

# get fits tdim shape for given array
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

# get fits tform string for given table type and data array.
fits_tform{T}(::Type{TableHDU}, A::Array{T}) = "$(prod(fits_tdim(A)))$(fits_tform_char(T))"

# For string arrays with 2+ dimensions, write tform as rAw. Otherwise,
# cfitsio doesn't recognize that multiple strings should be written to
# a single row, even if TDIM is set to 2+ dimensions.
fits_tform(::Type{TableHDU}, A::Vector{ASCIIString}) = "$(maximum(length, A))A"
fits_tform(::Type{TableHDU}, A::Array{ASCIIString}) = "$(prod(fits_tdim(A)))A$(maximum(length, A))"

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
    nrows_in = Array(Int64, 1)
    status = Cint[0]

    # fits_read_btblhdrll (Can pass NULL for return fields not needed.)
    ccall(("ffghbnll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Cint}, Ptr{Ptr{Uint8}},  # nrows, tfields, ttype
           Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}}, Ptr{Uint8},  # tform,tunit,extname
           Ptr{Clong}, Ptr{Cint}),  # pcount, status
          hdu.fitsfile.ptr, ncols, nrows_in, C_NULL, colnames_in, coltypes_in,
          C_NULL, C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    nrows = nrows_in[1]
    colnames = [bytestring(pointer(item)) for item in colnames_in]
    coltypes = [bytestring(pointer(item)) for item in coltypes_in]

    @printf io "file: %s\nextension: %d\ntype: BINARY TABLE\nrows: %d\ncolumns:" fits_file_name(hdu.fitsfile) hdu.ext nrows
    for i in 1:ncols
        @printf io "\n    %s (%s)" colnames[i] coltypes[i]
    end
end

function show(io::IO, hdu::ASCIITableHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    # allocate return arrays for column names & types
    colnames_in = [Array(Uint8, 70) for i=1:ncols]
    coltypes_in = [Array(Uint8, 70) for i=1:ncols]
    nrows_in = Array(Int64, 1)
    status = Cint[0]

    # fits_read_atblhdrll (Can pass NULL for return fields not needed)
    ccall(("ffghtbll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Int64}, Ptr{Cint},  # rowlen, nrows, tfields
           Ptr{Ptr{Uint8}}, Ptr{Clong}, Ptr{Ptr{Uint8}},  # ttype, tbcol, tform
           Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Cint}),  # tunit, extname, status
          hdu.fitsfile.ptr, ncols,
          C_NULL, nrows_in, C_NULL,
          colnames_in, C_NULL, coltypes_in,
          C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    nrows = nrows_in[1]
    colnames = [bytestring(pointer(item)) for item in colnames_in]
    coltypes = [bytestring(pointer(item)) for item in coltypes_in]

    @printf io "file: %s\nextension: %d\ntype: ASCII TABLE\nrows: %d\ncolumns:" fits_file_name(hdu.fitsfile) hdu.ext nrows
    for i in 1:ncols
        @printf io "\n    %s (%s)" colnames[i] coltypes[i]
    end
    print(io, "\n")
end

# Add a new TableHDU to a FITS object
function write_impl(f::FITS, colnames::Vector{ASCIIString}, coldata::Vector,
                    hdutype, extname, header, units)
    fits_assert_open(f.fitsfile)

    # move to last HDU; table will be added after the CHDU
    fits_movabs_hdu(f.fitsfile, @compat(Int(fits_get_num_hdus(f.fitsfile))))

    # create an array of tform strings (which we will create pointers to)
    tform_str = [fits_tform(hdutype, a) for a in coldata]

    ncols = length(coldata)
    ttype = [pointer(name) for name in colnames]
    tform = [pointer(s) for s in tform_str]

    # get units
    if isa(units, Nothing)
        tunit = C_NULL
    else
        tunit = Ptr{Uint8}[(haskey(units, name)? pointer(units[name]): C_NULL)
                           for name in colnames]
    end
    
    # extension name
    extname_ptr = (isa(extname, Nothing) ? convert(Ptr{Uint8}, C_NULL) :
                   pointer(extname))

    status = Cint[0]
    ccall(("ffcrtb", libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Cint, Ptr{Ptr{Uint8}}, Ptr{Ptr{Uint8}},
           Ptr{Ptr{Uint8}}, Ptr{Uint8}, Ptr{Cint}),
          f.fitsfile.ptr, table_type_code(hdutype), 0, ncols,  # 0 = nrows
          ttype, tform, tunit, extname_ptr, status)
    fits_assert_ok(status[1])

    # For binary tables, write tdim info
    if hdutype === TableHDU
        for (i, a) in enumerate(coldata)
            fits_write_tdim(f.fitsfile, i, fits_tdim(a))
        end
    end

    if isa(header, FITSHeader)
        write_header(f.fitsfile, header, true)
    end

    for (i, a) in enumerate(coldata)
        fits_write_col(f.fitsfile, i, 1, 1, a)
    end
    nothing
end

function write(f::FITS, colnames::Vector{ASCIIString}, coldata::Vector;
               units=nothing, header=nothing, hdutype=TableHDU,
               extname=nothing)
    if length(colnames) != length(coldata)
        error("length of colnames and length of coldata must match")
    end
    write_impl(f, colnames, coldata, hdutype, extname, header, units)
end

function write(f::FITS, data::Dict{ASCIIString};
               units=nothing, header=nothing, hdutype=TableHDU,
               extname=nothing)
    colnames = collect(keys(data))
    coldata = collect(values(data))
    write_impl(f, colnames, coldata, hdutype, extname, header, units)
end

# Read a table column into an array of the "equivalent type"
function read(hdu::Union(TableHDU, ASCIITableHDU), colname::ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname)

    # `eqcoltype`: do SCALE/ZERO conversion automatically
    typecode, repcnt, width = fits_get_eqcoltype(hdu.fitsfile, colnum)

    # BitArrays not yet supported.
    (abs(typecode) == 1) && error("BitArray ('X') columns not yet supported")

    T = CFITSIO_COLTYPE[typecode]

    # ASCII tables can only have scalar columns
    if isa(hdu, ASCIITableHDU)
        result = Array(T, nrows)
    else
        if abs(typecode) == 16
            # for strings, cfitsio only considers it to be a vector column if
            # width != repcnt, even if tdim is multi-valued.
            if repcnt == width
                result = Array(T, nrows)
            else
                rowsize = fits_read_tdim(hdu.fitsfile, colnum)
                # if rowsize isn't multi-valued, ignore it (we know it *is* a
                # vector column). If it is multi-valued, prefer it to repcnt,
                # width.
                if length(rowsize) == 1
                    result = Array(T, div(repcnt, width), nrows)
                else
                    result = Array(T, rowsize[2:end]..., nrows)
                end
            end
        else
            if repcnt == 1
                result = Array(T, nrows)
            else
                rowsize = fits_read_tdim(hdu.fitsfile, colnum)
                result = Array(T, rowsize..., nrows)
            end
        end
    end

    # TODO: allow altering first row and first element.
    fits_read_col(hdu.fitsfile, colnum, 1, 1, result)

    return result
end
