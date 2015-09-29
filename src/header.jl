# FITSHeader methods

# -----------------------------------------------------------------------------
# Helper functions
#
# Used here and in other files. Functions that operate on FITSFile
# start with `fits_`.

function try_parse_hdrval(::Type{Bool}, s::ASCIIString)
    if length(s) == 1
        if s[1] == 'T'
            return Nullable(true)
        elseif s[1] == 'F'
            return Nullable(false)
        end
    end
    return Nullable{Bool}()
end

# Note that trailing whitespace is not significant in FITS header
# keywords, but *leading* whitespace is, so "'    '" parses as " " (a
# single space).  See CFITSIO manual section 4.5 for details.
#
# TODO: parse '' within the string as a single '.
function try_parse_hdrval(::Type{ASCIIString}, s::ASCIIString)
    if length(s) < 2 || s[1] != '\'' || s[end] != '\''
        return Nullable{ASCIIString}()
    end

    i = endof(s) - 1
    while i > 2
        if s[i] != ' '
            return Nullable(s[2:i])
        end
        i -= 1
    end
    return Nullable(s[2:i])
end

try_parse_hdrval(::Type{Float64}, s::ASCIIString) = tryparse(Float64, s)

# hack for integers in Julia v0.3: tryparse(Int, s) not available in Compat.
if VERSION > v"0.4.0-dev+3864"
    try_parse_hdrval(::Type{Int}, s::ASCIIString) = tryparse(Int, s)
else
    try_parse_hdrval(::Type{Int}, s::ASCIIString) = try
        Nullable(parseint(s))
    catch e
        Nullable{Int}()
    end
end

# Try to parse the header value as any type
function try_parse_hdrval(s::ASCIIString)
    length(s) == 0 && return Nullable(nothing)

    nb = try_parse_hdrval(Bool, s)
    isnull(nb) || return nb

    ns = try_parse_hdrval(ASCIIString, s)
    isnull(ns) || return ns

    ni = try_parse_hdrval(Int, s)
    isnull(ni) || return ni

    nf = try_parse_hdrval(Float64, s)
    isnull(nf) || return nf

    return Nullable{Any}()
end

# functions for displaying header values in show(io, header)
hdrval_repr(v::Bool) = v ? "T" : "F"
hdrval_repr(v::ASCIIString) = @sprintf "'%s'" v
hdrval_repr(v::@compat(Union{AbstractFloat, Integer})) = string(v)

# returns one of: ASCIIString, Bool, Int, Float64, nothing
# (never error)
function parse_header_val(s::ASCIIString)
    nval = try_parse_hdrval(s)
    return isnull(nval) ? s : get(nval)
end

# Try to read the raw keys in order given; returns Nullable.
# (null if no key exists or if parsing an existing key is unsuccessful.)
function fits_try_read_keys{T}(f::FITSFile, ::Type{T}, keys)
    status = Cint[0]
    value = Array(@compat(UInt8), 71)
    for key in keys
        ccall((:ffgkey, libcfitsio), Cint,
              (Ptr{Void},Ptr{@compat(UInt8)},Ptr{@compat(UInt8)},Ptr{@compat(UInt8)},Ptr{Cint}),
              f.ptr, bytestring(key), value, C_NULL, status)

        # If the key is found, return it. If there was some other error
        # besides key not found, throw an error.
        if status[1] == 0
            return try_parse_hdrval(T, bytestring(pointer(value)))
        elseif status[1] != 202
            error(fits_get_errstatus(status[1]))
        end
    end
    return Nullable{T}()
end

# Build a string with extension keywords, if present.
# This is a helper function for show(::HDU).
const EXTNAME_KEYS = ["EXTNAME", "HDUNAME"]
const EXTVER_KEYS = ["EXTVER", "HDUVER"]
fits_try_read_extname(f::FITSFile) =
    fits_try_read_keys(f, ASCIIString, EXTNAME_KEYS)
fits_try_read_extver(f::FITSFile) = fits_try_read_keys(f, Int, EXTVER_KEYS)

function fits_get_ext_info_string(f::FITSFile)
    extname = fits_try_read_extname(f)
    extver = fits_try_read_extver(f)
    if !isnull(extname) && !isnull(extver)
        return " (name=$(repr(get(extname))), ver=$(get(extver)))"
    elseif !isnull(extname)
        return " (name=$(repr(get(extname))))"
    end
    return ""
end


# Return indicies of reserved keys in a header.
# This is more complex than you would think because some reserved keys
# are only reserved when other keys are present. Also, in general a key
# may appear more than once in a header.
const RESERVED_KEYS = ["SIMPLE","EXTEND","XTENSION","BITPIX","PCOUNT","GCOUNT",
                       "THEAP","EXTNAME","BUNIT","BSCALE","BZERO","BLANK",
                       "ZQUANTIZ","ZDITHER0","ZIMAGE","ZCMPTYPE","ZSIMPLE",
                       "ZTENSION","ZPCOUNT","ZGCOUNT","ZBITPIX","ZEXTEND",
                       "CHECKSUM","DATASUM"]
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


# Write header to CHDU.
# If `clean` is true, skip writing reserved header keywords.
function fits_write_header(f::FITSFile, hdr::FITSHeader, clean::Bool=true)
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
# Public API

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
    key = Array(@compat(UInt8), 81)
    value = Array(@compat(UInt8), 81)
    comment = Array(@compat(UInt8), 81)
    status = Cint[0]

    nkeys, morekeys = fits_get_hdrspace(hdu.fitsfile)

    # Initialize output arrays
    keys = Array(ASCIIString, nkeys)
    values = Array(Any, nkeys)
    comments = Array(ASCIIString, nkeys)
    for i=1:nkeys
        ccall((:ffgkyn,libcfitsio), Cint,
              (Ptr{Void},Cint,Ptr{@compat(UInt8)},Ptr{@compat(UInt8)},Ptr{@compat(UInt8)},Ptr{Cint}),
              hdu.fitsfile.ptr, i, key, value, comment, status)
        keys[i] = bytestring(pointer(key))
        values[i] = parse_header_val(bytestring(pointer(value)))
        comments[i] = bytestring(pointer(comment))
    end
    fits_assert_ok(status[1])
    FITSHeader(keys, values, comments)
end

function read_header(hdu::HDU, ::Type{ASCIIString})
    # Return the header as a raw string.

    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    fits_hdr2str(hdu.fitsfile)
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
                val = hdrval_repr(hdr.values[i])
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
