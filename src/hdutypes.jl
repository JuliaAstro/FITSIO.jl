
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

# -----------------------------------------------------------------------------
# FITS methods

function length(f::FITS)
    fits_assert_open(f.fitsfile)
    int(fits_get_num_hdus(f.fitsfile))
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
function _convert(val::String)
    try
        return int(val)
    catch
        try
            return float(val)
        catch
        end
    end
    val == "T" ? true :
    val == "F" ? false :
    val == "" ? nothing :
    length(val) > 1 && val[1] == '\'' && val[end] == '\'' ? val[2:end-1] :
    error("couldn't parse keyword value: \"$val\"")
end

function read_key(hdu::HDU, keynum::Integer)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    name, val, comment = fits_read_keyn(hdu.fitsfile, keynum)
    name, _convert(val), comment
end

function read_key(hdu::HDU, keyname::String)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    name, val = fits_read_keyword(hdu.fitsfile, heyname)
    name, _convert(val)
end

function read_header(hdu::HDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # Below, we use a direct call to ffgkyn so that we can keep reusing the
    # same keyname, value, comment buffers.
    keyname = Array(Uint8, 9)
    value = Array(Uint8, 71)
    comment = Array(Uint8, 71)
    status = Int32[0]

    nkeys, morekeys = fits_get_hdrspace(hdu.fitsfile)

    header = Array(Dict{String, Any}, nkeys)
    for i=1:nkeys
        ccall((:ffgkyn,libcfitsio), Int32,
              (Ptr{Void},Int32,Ptr{Uint8},Ptr{Uint8},Ptr{Uint8},Ptr{Int32}),
              hdu.fitsfile.ptr, i, keyname, value, comment, status)
        header[i] = ["name"=>bytestring(convert(Ptr{Uint8},keyname)),
                     "value"=>_convert(bytestring(convert(Ptr{Uint8},value))),
                     "comment"=>bytestring(convert(Ptr{Uint8},comment))]
    end
    fits_assert_ok(status[1])
    return header
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

# Read all or a subset of an HDU
function getindex(hdu::ImageHDU, rs::Range...)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # construct first, last and step vectors
    firsts = Clong[first(r) for r in rs]
    lasts = Clong[last(r) for r in rs]
    steps = Clong[step(r) for r in rs]

    # construct output array
    bitpix = fits_get_img_equivtype(hdu.fitsfile)
    
    datasz = [length(r) for r in rs]
    data = Array(bitpix_to_type[bitpix], datasz...)
    fits_read_subset(hdu.fitsfile, firsts, lasts, steps, data)
    data
end

# Add a new ImageHDU to a FITS object
# The following Julia data types are supported for writing images by cfitsio:
# Uint8, Int8, Uint16, Int16, Uint32, Int32, Int64, Float32, Float64
function write{T}(f::FITS, data::Array{T})
    fits_assert_open(f.fitsfile)
    s = size(data)
    fits_create_img(f.fitsfile, T, [s...])
    fits_write_pix(f.fitsfile, ones(Int, length(s)), length(data), data)
    nothing
end
