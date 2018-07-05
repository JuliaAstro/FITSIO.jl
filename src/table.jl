# TableHDU & ASCIITableHDU methods

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
                         (String,      'A',  16),
                         (UInt16,      'U',  20),
                         (Int16,       'I',  21),
                         (UInt32,      'V',  40),
                         (Int32,       'J',  41),
                         (Int64,       'K',  81),
                         (Float32,     'E',  42),
                         (Float64,     'D',  82),
                         (ComplexF32,  'C',  83),
                         (ComplexF64,  'M', 163))
    @eval fits_tform_char(::Type{$T}) = $tform
    CFITSIO_COLTYPE[code] = T
end
const FITSTableScalar = Union{UInt8, Int8, Bool, UInt16, Int16, UInt32,
                              Int32, Int64, Float32, Float64, ComplexF32,
                              ComplexF64}

# Helper function for reading information about a (binary) table column
# Returns: (eltype, rowsize, isvariable)
function fits_get_col_info(f::FITSFile, colnum::Integer)
    eqtypecode, repeat, width = fits_get_eqcoltype(f, colnum)
    isvariable = eqtypecode < 0
    eqtypecode = abs(eqtypecode)

    (eqtypecode == 1) && error("BitArray ('X') columns not yet supported")

    T = CFITSIO_COLTYPE[eqtypecode]

    if isvariable
        if T !== String
            T = Vector{T}
        end
        rowsize = Int[]
    else
        if T === String
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

# Parse max length from tform for a variable column
function var_col_maxlen(tform::String)
    maxlen = -1
    i = search(tform, '(')
    if i > 0
        j = search(tform, ')', i)
        if j > 0
            try maxlen = parseint(tform[i+1:j-1]) catch; end
        end
    end
    return maxlen
end

# Helper function for getting fits tdim shape for given array
fits_tdim(A::Array) = (ndims(A) == 1) ? [1] : [size(A, i) for i=1:ndims(A)-1]
function fits_tdim(A::Array{String})
    n = ndims(A)
    tdim = Vector{Int}(n)
    tdim[1] = maximum(length, A)
    for i=2:n
        tdim[n] = size(A, n-1)
    end
    tdim
end

# Helper function for getting fits tform string for given table type
# and data array.
fits_tform(::Type{TableHDU}, A::Array{T}) where {T} = "$(prod(fits_tdim(A)))$(fits_tform_char(T))"

# For string arrays with 2+ dimensions, write tform as rAw. Otherwise,
# cfitsio doesn't recognize that multiple strings should be written to
# a single row, even if TDIM is set to 2+ dimensions.
fits_tform(::Type{TableHDU}, A::Vector{String}) = "$(maximum(length, A))A"
fits_tform(::Type{TableHDU}, A::Array{String}) = "$(prod(fits_tdim(A)))A$(maximum(length, A))"

# variable length columns
fits_tform_v(::Type{TableHDU}, A::Vector{Vector{T}}) where {T<:FITSTableScalar} = "1P$(fits_tform_char(T))($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{String}) = "1PA($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{Vector}) = error("column data must be a leaf type: e.g., Vector{Vector{Int}}, not Vector{Vector{T}}.")
fits_tform_v(::Type{TableHDU}, ::Any) = error("variable length columns only supported for arrays of arrays and arrays of String")
fits_tform_v(::Type{ASCIITableHDU}, ::Any) = error("variable length columns not supported in ASCII tables")


fits_tform(::Type{ASCIITableHDU}, ::Vector{Int16}) = "I7"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Int32}) = "I12"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float32}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float64}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, A::Vector{String}) = "A$(maximum(length, A))"
fits_tform(::Type{ASCIITableHDU}, A::Vector) = error("unsupported type: $(eltype(A))")
fits_tform(::Type{ASCIITableHDU}, A::Array) = error("only 1-d arrays supported: dimensions are $(size(A))")

# for passing to fits_create_tbl.
table_type_code(::Type{ASCIITableHDU}) = Cint(1)
table_type_code(::Type{TableHDU}) = Cint(2)

function fits_read_table_header!(hdu::TableHDU, ncols, nrows,
                                 colnames_in, coltforms_in, status)
    # fits_read_btblhdrll (Can pass NULL for return fields not needed.)
    ccall(("ffghbnll", libcfitsio), Cint,
          (Ptr{Cvoid}, Cint,  # Inputs: fitsfile, maxdim
           Ref{Int64}, Ptr{Cint}, Ptr{Ptr{UInt8}},  # nrows, tfields, ttype
           Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{UInt8},  # tform,tunit,extname
           Ptr{Clong}, Ref{Cint}),  # pcount, status
          hdu.fitsfile.ptr, ncols, nrows, C_NULL, colnames_in, coltforms_in,
          C_NULL, C_NULL, C_NULL, status)
end

function fits_read_table_header!(hdu::ASCIITableHDU, ncols, nrows,
                                 colnames_in, coltforms_in, status)
    # fits_read_atblhdrll (Can pass NULL for return fields not needed)
    ccall(("ffghtbll", libcfitsio), Cint,
          (Ptr{Cvoid}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ref{Int64}, Ptr{Cint},  # rowlen, nrows, tfields
           Ptr{Ptr{UInt8}}, Ptr{Clong}, Ptr{Ptr{UInt8}},  # ttype, tbcol, tform
           Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ref{Cint}),  # tunit, extname, status
          hdu.fitsfile.ptr, ncols, C_NULL, nrows, C_NULL, colnames_in,
          C_NULL, coltforms_in, C_NULL, C_NULL, status)
end

function columns_names_tforms(hdu::Union{ASCIITableHDU,TableHDU})
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    # allocate return arrays for column names & types
    colnames_in  = [Vector{UInt8}(undef, 70) for i=1:ncols]
    coltforms_in = [Vector{UInt8}(undef, 70) for i=1:ncols]
    nrows = Ref{Int64}()
    status = Ref{Cint}(0)

    fits_read_table_header!(hdu, ncols, nrows, colnames_in, coltforms_in, status)
    fits_assert_ok(status[])

    # parse out results
    colnames = [unsafe_string(pointer(item)) for item in colnames_in]
    coltforms = [unsafe_string(pointer(item)) for item in coltforms_in]
    return colnames, coltforms, ncols, nrows
end

"""
    colnames(hdu) -> Vector{String}

Return the names of columns in a table HDU.
"""
colnames(hdu::Union{ASCIITableHDU,TableHDU}) = columns_names_tforms(hdu)[1]

function show(io::IO, hdu::TableHDU)
    colnames, coltforms, ncols, nrows = columns_names_tforms(hdu)
    # get some more information for all the columns
    coltypes    = Vector{String}(ncols)
    colrowsizes = Vector{String}(ncols)
    showlegend = false
    for i in 1:ncols
        T, rowsize, isvariable = fits_get_col_info(hdu.fitsfile, i)
        coltypes[i] = repr(T)
        colrowsizes[i] = length(rowsize) > 0 ? repr(tuple(rowsize...)) : ""
        if isvariable
            coltypes[i] *= "*"
            showlegend = true
        end
    end

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(hdu.ext)$(fits_get_ext_info_string(hdu.fitsfile))
    Type: Table
    Rows: $(nrows[])
    Columns: """)
    show_ascii_table(
        io, ["Name", "Size", "Type", "TFORM"],
        Vector{String}[colnames, colrowsizes, coltypes, coltforms],
        2, 9)
    if showlegend
        print(io, "\n\n         (*) = variable-length column")
    end
end

function show(io::IO, hdu::ASCIITableHDU)
    colnames, coltforms, ncols, nrows = columns_names_tforms(hdu)
    # Get additional info
    coltypes = Vector{String}(ncols)
    for i in 1:ncols
        eqtypecode, repeat, width = fits_get_eqcoltype(hdu.fitsfile, i)
        T = CFITSIO_COLTYPE[eqtypecode]
        coltypes[i] = repr(T)
    end

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(hdu.ext)$(fits_get_ext_info_string(hdu.fitsfile))
    Type: ASCIITable
    Rows: $(nrows[])
    Columns: """)
    show_ascii_table(io, ["Name", "Type", "TFORM"],
                     Vector{String}[colnames, coltypes, coltforms], 2, 9)
end

# Write a variable length array column of numbers
# (separate implementation from normal fits_write_col function because
#  we must make separate calls to `fits_write_col` for each row.)
function fits_write_var_col(f::FITSFile, colnum::Integer,
                            data::Vector{Vector{T}}) where T
    for i=1:length(data)
        fits_write_col(f, colnum, i, 1, data[i])
    end
end

function fits_write_var_col(f::FITSFile, colnum::Integer,
                            data::Vector{String})
    for el in data; fits_assert_isascii(el); end
    status = Ref{Cint}(0)
    buffer = Ref{Ptr{UInt8}}()  # holds the address of the current row
    for i=1:length(data)
        buffer[] = pointer(data[i])

        # Note that when writing to a variable ASCII column, the
        # ‘firstelem’ and ‘nelements’ parameter values in the
        # fits_write_col routine are ignored and the number of
        # characters to write is simply determined by the length of
        # the input null-terminated character string.
        ccall((:ffpcls, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Int64, Int64, Int64, Ref{Ptr{UInt8}},
               Ref{Cint}),
              f.ptr, colnum, i, 1, length(data[i]), buffer, status)
        fits_assert_ok(status[])
    end
end

# Add a new TableHDU to a FITS object
function write_internal(f::FITS, colnames::Vector{String},
                        coldata::Vector, hdutype, name, ver, header, units,
                        varcols)
    fits_assert_open(f.fitsfile)
    for el in colnames; fits_assert_isascii(el); end

    # move to last HDU; table will be added after the CHDU
    nhdus = Int(fits_get_num_hdus(f.fitsfile))
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
    tform_str = Vector{String}(ncols)
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
        tunit = Ptr{UInt8}[(haskey(units, n) ? pointer(units[n]) : C_NULL)
                           for n in colnames]
    end

    # extension name
    name_ptr = (isa(name, Nothing) ? Ptr{UInt8}(C_NULL) :
                   pointer(name))

    status = Ref{Cint}(0)
    ccall(("ffcrtb", libcfitsio), Cint,
          (Ptr{Cvoid}, Cint, Int64, Cint, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}},
           Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ref{Cint}),
          f.fitsfile.ptr, table_type_code(hdutype), 0, ncols,  # 0 = nrows
          ttype, tform, tunit, name_ptr, status)
    fits_assert_ok(status[])

    # For binary tables, write tdim info
    if hdutype === TableHDU
        for (i, a) in enumerate(coldata)
            isvarcol[i] || fits_write_tdim(f.fitsfile, i, fits_tdim(a))
        end
    end

    if isa(header, FITSHeader)
        fits_write_header(f.fitsfile, header, true)
    end
    if isa(ver, Integer)
        fits_update_key(f.fitsfile, "EXTVER", ver)
    end

    for (i, a) in enumerate(coldata)
        if isvarcol[i]
            fits_write_var_col(f.fitsfile, i, a)
        else
            fits_write_col(f.fitsfile, i, 1, 1, a)
        end
    end
    nothing
end


"""
    write(f::FITS, colnames, coldata; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)

Same as `write(f::FITS, data::Dict; ...)` but providing column names
and column data as a separate arrays. This is useful for specifying
the order of the columns. Column names must be `Vector{String}`
and column data must be a vector of arrays.
"""
function write(f::FITS, colnames::Vector{String}, coldata::Vector;
               units=nothing, header=nothing, hdutype=TableHDU,
               name=nothing, ver=nothing, varcols=nothing)
    if length(colnames) != length(coldata)
        error("length of colnames and length of coldata must match")
    end
    write_internal(f, colnames, coldata, hdutype, name, ver, header,
                   units, varcols)
end


"""
    write(f::FITS, data::Dict; hdutype=TableHDU, name=nothing, ver=nothing, header=nothing, units=nothing, varcols=nothing)

Create a new table extension and write data to it. If the FITS file is
currently empty then a dummy primary array will be created before
appending the table extension to it. `data` should be a dictionary
with String keys (giving the column names) and Array values
(giving data to write to each column). The following types are
supported in binary tables: `UInt8`, `Int8`, `UInt16`, `Int16`,
`UInt32`, `Int32`, `Int64`, `Float32`, `Float64`, `Complex{Float32}`,
`Complex{Float64}`, `String`, `Bool`.

Optional inputs:

- `hdutype`: Type of table extension to create. Can be either
  `TableHDU` (binary table) or `ASCIITableHDU` (ASCII table).
- `name`: Name of extension.
- `ver`: Version of extension (Int).
- `header`: FITSHeader instance to write to new extension.
- `units`: Dictionary mapping column name to units (as a string).
- `varcols`: An array giving the column names or column indicies to
  write as "variable-length columns".

!!! note "Variable length columns"

    Variable length columns allow a column's row entries to contain
    arrays of different lengths. They can potentially save diskspace
    when the rows of a column vary greatly in length, as the column
    data is all written to a contiguous heap area at the end of the
    table. Only column data of type `Vector{String}` or types
    such as `Vector{Vector{UInt8}}` can be written as variable
    length columns. In the second case, ensure that the column data
    type is a *leaf type*. That is, the type cannot be
    `Vector{Vector{T}}`, which would be an array of arrays having
    potentially non-uniform element types (which would not be writable
    as a FITS table column).
"""
function write(f::FITS, data::Dict{String};
               units=nothing, header=nothing, hdutype=TableHDU,
               name=nothing, ver=nothing, varcols=nothing)
    colnames = collect(keys(data))
    coldata = collect(values(data))
    write_internal(f, colnames, coldata, hdutype, name, ver, header,
                   units, varcols)
end

# Read a variable length array column of numbers
# (separate implementation from normal fits_read_col function because
# the length of each vector must be determined for each row.
function fits_read_var_col(f::FITSFile, colnum::Integer,
                           data::Vector{Vector{T}}) where T
    nrows = length(data)
    for i=1:nrows
        repeat, offset = fits_read_descript(f, colnum, i)
        data[i] = Vector{T}(repeat)
        fits_read_col(f, colnum, i, 1, data[i])
    end
end

# Read a variable length array column of strings
# (Must be separate implementation from normal fits_read_col function because
# the length of each string must be determined for each row.)
function fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{String})
    status = Ref{Cint}(0)
    bufptr = Ref{Ptr{UInt8}}()  # holds a pointer to the current row buffer
    for i=1:length(data)
        repeat, offset = fits_read_descript(f, colnum, i)
        buffer = Vector{UInt8}(repeat)
        bufptr[] = pointer(buffer)
        ccall((:ffgcvs, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Int64, Int64, Int64,
               Ptr{UInt8}, Ref{Ptr{UInt8}}, Ptr{Cint}, Ref{Cint}),
              f.ptr, colnum, i, 1, repeat, " ", bufptr, C_NULL, status)
        fits_assert_ok(status[])

        # Create string out of the buffer, terminating at null characters
        zeropos = search(buffer, 0x00)
        data[i] = (zeropos >= 1) ? String(buffer[1:(zeropos-1)]) :
                                   String(buffer)
    end
end

# Read a table column
function read(hdu::ASCIITableHDU, colname::String; case_sensitive::Bool=true)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=case_sensitive)

    typecode, repcnt, width = fits_get_eqcoltype(hdu.fitsfile, colnum)
    T = CFITSIO_COLTYPE[typecode]

    result = Vector{T}(nrows)
    fits_read_col(hdu.fitsfile, colnum, 1, 1, result)

    return result
end


"""
    read(hdu, colname; case_sensitive=true)

Read a column as an array from the given table HDU.

The column name may contain wild card characters (`*`, `?`, or
`#`). The `*` wild card character matches any sequence of
characters (including zero characters) and the `?` character
matches any single character. The `#` wildcard will match any
consecutive string of decimal digits (0-9). The string must match a
unique column.  The optional boolean keyword `case_sensitive`,
`true` by default, specifies whether the column name is to be
considered case sensitive.
"""
function read(hdu::TableHDU, colname::String; case_sensitive::Bool=true)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=case_sensitive)

    T, rowsize, isvariable = fits_get_col_info(hdu.fitsfile, colnum)

    result = Array{T}(rowsize..., nrows)

    if isvariable
        fits_read_var_col(hdu.fitsfile, colnum, result)
    else
        fits_read_col(hdu.fitsfile, colnum, 1, 1, result)
    end

    return result
end
