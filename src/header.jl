# FITSHeader methods

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
