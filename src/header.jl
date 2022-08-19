# FITSHeader methods

# -----------------------------------------------------------------------------
# Helper functions
#
# Used here and in other files. Functions that operate on FITSFile
# start with `fits_`.

"""
    try_parse_hdrval(T, s) -> x::Union{T,Nothing}

Attempt to parse string `s` for a FITS card value of type `T` (`String`,
`Bool`, `Int` or `Float64`) and yield either a value of type `T` or
`nothing` if parsing is unsuccessful.

See also: [`parse_header_val`](@ref).

"""
function try_parse_hdrval(::Type{Bool}, s::String)
    if length(s) == 1
        if s[1] == 'T'
            return true
        elseif s[1] == 'F'
            return false
        end
    end
    return nothing
end

# Note that trailing whitespaces are not significant in FITS header
# keywords, but *leading* whitespaces are, so "'    '" parses as " " (a
# single space).  See CFITSIO manual section 4.5 for details.
#
# TODO: parse '' within the string as a single '.
function try_parse_hdrval(::Type{String}, s::String)
    if length(s) < 2 || s[1] != '\'' || s[end] != '\''
        return nothing
    end

    i = lastindex(s) - 1
    while i > 2
        if s[i] != ' '
            return s[2:i]
        end
        i -= 1
    end
    return s[2:i]
end

try_parse_hdrval(::Type{T}, s::String) where {T<:Union{Int,Float64}} =
    tryparse(T, s)

# functions for displaying header values in show(io, header)
hdrval_repr(v::Bool) = v ? "T" : "F"
hdrval_repr(v::String) = @sprintf "'%-8s'" v
hdrval_repr(v::Union{AbstractFloat, Integer}) = string(v)

"""
    parse_header_val(s::String) -> x::Union{String,Bool,Int,Float64,Nothing}

Parse the FITS card value in the string `s` and return a value of type
`String`, `Bool`, `Int`, `Float64` or `Nothing`.  The latter indicates an
empty value.  This method never throws an error; if the value cannot be
parsed it is returned as it.

See also: [`try_parse_hdrval`](@ref).
"""
function parse_header_val(s::String)
    length(s) == 0 && return nothing

    vb = try_parse_hdrval(Bool, s)
    vb === nothing || return vb

    vs = try_parse_hdrval(String, s)
    vs === nothing || return vs

    vi = try_parse_hdrval(Int, s)
    vi === nothing || return vi

    vf = try_parse_hdrval(Float64, s)
    vf === nothing || return vf

    # FIXME: error("invalid FITS keyword value")
    return s
end

"""
    fits_try_read_keys(f::FITSFile, ::Type{T}, keys)

Try to read the raw FITS keys in given order if FITS handle `f` and
return a value of type `T` or [`nothing`](@ref) if no key exists or if
parsing an existing key is unsuccessful.

See also: [`try_parse_hdrval`](@ref).

"""
function fits_try_read_keys(f::FITSFile, ::Type{T}, keys) where T
    status = Cint[0]
    value = Vector{UInt8}(undef, 71)
    for key in keys
        ccall((:ffgkey, libcfitsio), Cint,
              (Ptr{Cvoid},Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint}),
              f.ptr, key, value, C_NULL, status)

        # If the key is found, return it. If there was some other error
        # besides key not found, throw an error.
        if status[1] == 0
            return try_parse_hdrval(T, unsafe_string(pointer(value)))
        elseif status[1] != 202
            error(fits_get_errstatus(status[1]))
        end
    end
    return nothing
end

# Build a string with extension keywords, if present.
# This is a helper function for show(::HDU).
const EXTNAME_KEYS = ["EXTNAME", "HDUNAME"]
const EXTVER_KEYS = ["EXTVER", "HDUVER"]
fits_try_read_extname(f::FITSFile) :: Union{String,Nothing} =
    fits_try_read_keys(f, String, EXTNAME_KEYS)
fits_try_read_extver(f::FITSFile) :: Union{Int,Nothing} =
    fits_try_read_keys(f, Int, EXTVER_KEYS)

function fits_get_ext_info_string(f::FITSFile)
    extname = fits_try_read_extname(f)
    extver = fits_try_read_extver(f)
    if extname === nothing
        return ""
    elseif extver === nothing
        return " (name=$extname)"
    else
        return " (name=$extname, ver=$extver)"
    end
end

# Return indices of reserved keys in a header.
# This is more complex than you would think because some reserved keys
# are only reserved when other keys are present. Also, in general a key
# may appear more than once in a header.
const RESERVED_KEYS = ["SIMPLE","EXTEND","XTENSION","BITPIX","PCOUNT","GCOUNT",
                       "THEAP","EXTNAME","BSCALE","BZERO","BLANK",
                       "ZQUANTIZ","ZDITHER0","ZIMAGE","ZCMPTYPE","ZSIMPLE",
                       "ZTENSION","ZPCOUNT","ZGCOUNT","ZBITPIX","ZEXTEND",
                       "CHECKSUM","DATASUM"]
function reserved_key_indices(hdr::FITSHeader)
    nhdr = length(hdr)
    indices = Int[]
    for i=1:nhdr
        if in(hdr.keys[i], RESERVED_KEYS)
            push!(indices, i)
        end
    end

    # Note that this removes anything matching NAXIS\d regardless of # of axes.
    if in("NAXIS", hdr.keys)
        for i=1:nhdr
            if occursin(r"^NAXIS\d*$", hdr.keys[i])
                push!(indices, i)
            end
        end
    end

    if in("ZNAXIS", hdr.keys)
        for i=1:nhdr
            if (occursin(r"^ZNAXIS\d*$", hdr.keys[i]) ||
                occursin(r"^ZTILE\d*$", hdr.keys[i]) ||
                occursin(r"^ZNAME\d*$", hdr.keys[i]) ||
                occursin(r"^ZVAL\d*$", hdr.keys[i]))
                push!(indices, i)
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
                if occursin(re, hdr.keys[i])
                    push!(indices, i)
                end
            end
        end
    end

    return indices
end


# Write header to CHDU.
# If `clean` is true, skip writing reserved header keywords.
function fits_write_header(f::FITSFile, hdr::FITSHeader, clean::Bool=true)
    indices = clean ? reserved_key_indices(hdr) : Int[]
    for i=1:length(hdr)
        if clean && in(i, indices)
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

"""
    read_key(hdu::HDU, key::String) -> (value, comment)

Read the HDU header record specified by keyword and return a tuple where
`value` is the keyword parsed value (of type `String`, `Bool`, `Int`,
`Float64` or `Nothing`), `comment` is the keyword comment (as a string).
Throw an error if `key` is not found.

    read_key(hdu::HDU, key::Integer) -> (keyname, value, comment)

Same as above but FITS card is specified by its position and returns a 3
element tuple where `keyname` is the keyword name (a string).
"""
function read_key(hdu::HDU, key::Integer)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    keyout, value, comment = fits_read_keyn(hdu.fitsfile, key)
    keyout, parse_header_val(value), comment
end

function read_key(hdu::HDU, key::String)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    value, comment = fits_read_keyword(hdu.fitsfile, key)
    parse_header_val(value), comment
end


"""
    write_key(hdu::HDU, key::String, value[, comment])

Write a keyword value the HDU's header. `value` can be a standard
header type (`String`, `Bool`, `Integer`, `AbstractFloat`) or
`nothing`, in which case the value part of the record will be
empty. If the keyword already exists, the value will be
overwritten. The comment will only be overwritten if given. If the
keyword does not already exist, a new record will be appended at the
end of the header.
"""
function write_key(hdu::HDU, key::String,
                   value::Union{String, Bool, Integer, AbstractFloat, Nothing},
                   comment::Union{String, Ptr{Cvoid}}=C_NULL)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_update_key(hdu.fitsfile, key, value, comment)
end





"""
    read_header(filename::AbstractString, hduindex = 1) -> FITSHeader

Convenience function to read the entire header corresponding to the HDU at index `hduindex` contained
in the FITS file named `filename`. Functionally `read_header(filename, hduindex)` is equivalent to

```julia
FITS(filename, "r") do f
    read_header(f[hduindex])
end
```
"""
function read_header(filename::AbstractString, hduindex = 1)
    FITS(filename, "r") do f
        read_header(f[hduindex])
    end
end

"""
    read_header(hdu::HDU) -> FITSHeader

Read the entire header from the given HDU and return a `FITSHeader` object.
The value of each header record is parsed as `Int`, `Float64`, `String`,
`Bool` or `nothing` according to the FITS standard.

If the value cannot be parsed according to the FITS standard, the value is
stored as the raw unparsed `String`.
"""
function read_header(hdu::HDU)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    # Below, we use a direct call to ffgkyn so that we can keep reusing the
    # same buffers.
    key = Vector{UInt8}(undef, 81)
    value = Vector{UInt8}(undef, 81)
    comment = Vector{UInt8}(undef, 81)
    status = Cint[0]

    nkeys, morekeys = fits_get_hdrspace(hdu.fitsfile)

    # Initialize output arrays
    keys = Vector{String}(undef, nkeys)
    values = Vector{Any}(undef, nkeys)
    comments = Vector{String}(undef, nkeys)
    for i=1:nkeys
        ccall((:ffgkyn,libcfitsio), Cint,
              (Ptr{Cvoid},Cint,Ptr{UInt8},Ptr{UInt8},Ptr{UInt8},Ptr{Cint}),
              hdu.fitsfile.ptr, i, key, value, comment, status)
        keys[i] = unsafe_string(pointer(key))
        values[i] = parse_header_val(unsafe_string(pointer(value)))
        comments[i] = unsafe_string(pointer(comment))
    end
    fits_assert_ok(status[1])
    FITSHeader(keys, values, comments)
end


"""
    read_header(hdu::HDU, String) -> String

Read the entire header from the given HDU as a single string.
"""
function read_header(hdu::HDU, ::Type{String})
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    fits_hdr2str(hdu.fitsfile)
end


"""
    length(hdr::FITSHeader)

Number of records in header of HDU.
"""
length(hdr::FITSHeader) = length(hdr.keys)

"""
    haskey(hdr::FITSHeader, key::String)

Returns true if `key` exists in header, otherwise false.
"""
haskey(hdr::FITSHeader, key::String) = in(key, hdr.keys)

"""
    keys(hdr::FITSHeader)

Array of keywords in header of HDU (not a copy).
"""
keys(hdr::FITSHeader) = hdr.keys

"""
    values(hdr::FITSHeader)

Array of values in header of HDU (not a copy).
"""
values(hdr::FITSHeader) = hdr.values

getkey(hdr::FITSHeader, key::String, default) = haskey(hdr, key) ? key : default
get(hdr::FITSHeader, key::String, default) = haskey(hdr, key) ? hdr[key] : default
get(f::Function, hdr::FITSHeader, key::String) = haskey(hdr, key) ? hdr[key] : f()

getindex(hdr::FITSHeader, key::String) = hdr.values[hdr.map[key]]
getindex(hdr::FITSHeader, i::Integer) = hdr.values[i]

function setindex!(hdr::FITSHeader, value::Any, key::String)
    fits_assert_isascii(key)
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


"""
    delete!(hdr::FITSHeader, key::String)

Delete a key in a FITS header. 
"""
function delete!(hdr::FITSHeader, key::String)

    # Throw error on missing key
    haskey(hdr, key) || throw(KeyError(key))

    index = hdr.map[key]

    # Delete the entries
    deleteat!(hdr.values, index)
    deleteat!(hdr.keys, index)
    deleteat!(hdr.comments, index)


    # Then re-create the map array
    map = Dict(zip(hdr.keys, 1:length(hdr.keys)))

    return nothing
end



# Comments
"""
    get_comment(hdr::FITSHeader, key_or_index::Union{String,Integer})

Get the comment based on keyword or index.
"""
get_comment(hdr::FITSHeader, key::String) = hdr.comments[hdr.map[key]]
get_comment(hdr::FITSHeader, i::Integer) = hdr.comments[i]

"""
    set_comment!(hdr::FITSHeader, key_or_index::Union{String,Integer}, comment::String)

Set the comment based on keyword or index.
"""
function set_comment!(hdr::FITSHeader, key::String, comment::String)
    fits_assert_isascii(comment)
    hdr.comments[hdr.map[key]] = comment
end
function set_comment!(hdr::FITSHeader, i::Integer, comment::String)
    fits_assert_isascii(comment)
    hdr.comments[i] = comment
end

# Display the header
function show(io::IO, hdr::FITSHeader)
    n = length(hdr)
    for i=1:n
        if hdr.keys[i] == "COMMENT" || hdr.keys[i] == "HISTORY"
            @printf io "%s %s" hdr.keys[i] hdr.comments[i][1:min(71, end)]
        else
            @printf io "%-8s" hdr.keys[i]
            if hdr.values[i] === nothing
                print(io, "                      ")
                rc = 50  # remaining characters on line
            elseif hdr.values[i] isa String
                val = hdrval_repr(hdr.values[i])
                @printf io "= %-20s" val
                rc = length(val) <= 20 ? 50 : 70 - length(val)
            else
                val = hdrval_repr(hdr.values[i])
                @printf io "= %20s" val
                rc = length(val) <= 20 ? 50 : 70 - length(val)
            end

            if length(hdr.comments[i]) > 0
                @printf io " / %s" hdr.comments[i][1:min(rc-3, end)]
            end
        end
        i != n && println(io)
    end
end

"""
    default_header(data::AbstractArray)

Creates a default header for the given array with the `SIMPLE`, `BITPIX`, `NAXIS`, `NAXIS*`, and `EXTEND` entries.
"""
function default_header(data::AbstractArray{T}) where T <: Number
    # assigning keys
    hdu_keys = ["SIMPLE",
                "BITPIX",
                "NAXIS",
                ("NAXIS$i" for i in 1:ndims(data))...,
                "EXTEND"]

    # assiging values
    hdu_values = [true,                                           # SIMPLE
                  CFITSIO.bitpix_from_type(T),                    # BITPIX
                  ndims(data),                                    # NAXIS
                  reverse(size(data))...,                         # size of each axis
                  true]                                           # EXTEND

    # assigning comments
    comments = ["file does conform to FITS standard",                                   # comment for SIMPLE
                "number of bits per data pixel",                                        # comment for BITPIX
                "number of data axes",                                                  # comment for NAXIS
                ("length of data axis $i" for i in 1:ndims(data))...,                   # comments for axis length
                "FITS dataset may contain extensions"]                                  # comment for EXTEND

    return FITSHeader(hdu_keys, hdu_values, comments)
end
