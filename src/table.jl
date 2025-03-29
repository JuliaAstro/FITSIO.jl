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
                         (UInt64,      'W',  80),
                         (Int64,       'K',  81),
                         (Float32,     'E',  42),
                         (Float64,     'D',  82),
                         (ComplexF32,  'C',  83),
                         (ComplexF64,  'M', 163))
    @eval fits_tform_char(::Type{$T}) = $tform
    CFITSIO_COLTYPE[code] = T
end
const FITSTableScalar = Union{UInt8, Int8, Bool, UInt16, Int16, UInt32,
                              Int32, UInt64, Int64, Float32, Float64, ComplexF32,
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

function safe_table_string_convert(buffer::Vector{UInt8})
    GC.@preserve buffer begin
        i = findfirst(iszero, buffer)
        isnothing(i) && return String(buffer)
        return String(view(buffer, firstindex(buffer):i-1))
    end
end
# Parse max length from tform for a variable column
function var_col_maxlen(tform::String)
    maxlen = -1
    i = something(findfirst(isequal('('), tform), 0)
    if i > 0
        j = something(findnext(isequal(')'), tform, i), 0)
        if j > i
            # FIXME: use `tryparse`
            try maxlen = parseint(tform[i+1:j-1]) catch end
        end
    end
    return maxlen
end

# Helper function for getting fits tdim shape for given array
fits_tdim(A::Array) = (ndims(A) == 1) ? [1] : [size(A, i) for i=1:ndims(A)-1]
function fits_tdim(A::Array{String})
    n = ndims(A)
    tdim = Vector{Int}(undef, n)
    tdim[1] = maximum(length, A)
    for i=2:n
        tdim[n] = size(A, n-1)
    end
    tdim
end

# Helper function for getting fits tform string for given table type
# and data array.
fits_tform(::Type{TableHDU}, A::Array{T}) where {T} =
    "$(prod(fits_tdim(A)))$(fits_tform_char(T))"

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
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    GC.@preserve hdu begin
        # Allocate return arrays for column names & types
        colnames_in  = [Vector{UInt8}(undef, 70) for i=1:ncols]
        coltforms_in = [Vector{UInt8}(undef, 70) for i=1:ncols]
        nrows = Ref{Int64}()
        status = Ref{Cint}(0)

        fits_read_table_header!(hdu, ncols, nrows, colnames_in, coltforms_in, status)
        fits_assert_ok(status[])

        # Parse out results safely
        colnames = [safe_table_string_convert(item) for item in colnames_in]
        coltforms = [safe_table_string_convert(item) for item in coltforms_in]
        return colnames, coltforms, ncols, nrows
    end
end

"""
    colnames(hdu) -> Vector{String}

Return the names of columns in a table HDU.
"""
colnames(hdu::Union{ASCIITableHDU,TableHDU}) = columns_names_tforms(hdu)[1]

function show(io::IO, hdu::TableHDU)
    if isdeleted(hdu)
        print(io, "Deleted HDU")
        return
    end
    colnames, coltforms, ncols, nrows = columns_names_tforms(hdu)
    # get some more information for all the columns
    coltypes    = Vector{String}(undef, ncols)
    colrowsizes = Vector{String}(undef, ncols)
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
    if isdeleted(hdu)
        print(io, "Deleted HDU")
        return
    end
    colnames, coltforms, ncols, nrows = columns_names_tforms(hdu)
    # Get additional info
    coltypes = Vector{String}(undef, ncols)
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
"""
    VarColHandler

Abstract type for handling variable length column operations in FITS files.
"""
abstract type VarColHandler end


"""
    StringVarColHandler <: VarColHandler

Handler for variable length string columns in FITS files.
Used to handle reading and writing of variable-length string arrays.
"""
struct StringVarColHandler <: VarColHandler end

"""
    NumericVarColHandler <: VarColHandler

Handler for variable length numeric columns in FITS files.
Used to handle reading and writing of variable-length numeric arrays.
"""
struct NumericVarColHandler <: VarColHandler end
"""
    UnsupportedVarColHandler <: VarColHandler

Handler for unsupported variable length column types in FITS files.
Used to handle error cases for unsupported data types.
"""
struct UnsupportedVarColHandler <: VarColHandler end

"""
    varcolhandler(::Type{String})
    varcolhandler(::Type{Vector{T}}) where T <: Number
    varcolhandler(::Type{T}) where T
    varcolhandler(::Type{<:Vector{String}})
    varcolhandler(::Type{<:Vector{Vector{T}}} where T)

Determine the appropriate variable column handler for a given type.

# Returns
- `StringVarColHandler` for String types
- `NumericVarColHandler` for numeric vector types
- `UnsupportedVarColHandler` for unsupported types

# Examples
```julia
varcolhandler(String)                  # Returns StringVarColHandler()
varcolhandler(Vector{Float64})         # Returns NumericVarColHandler()
varcolhandler(Vector{String})          # Returns StringVarColHandler()
varcolhandler(Vector{Vector{Float64}}) # Returns NumericVarColHandler()
varcolhandler(Complex{Float64})        # Returns UnsupportedVarColHandler()
```
"""
function varcolhandler(::Type{String})
    StringVarColHandler()
end

function varcolhandler(::Type{Vector{T}}) where T <: Number
    NumericVarColHandler()
end

function varcolhandler(::Type{T}) where T
    UnsupportedVarColHandler()
end

# Type trait definitions
varcolhandler(::Type{<:Vector{String}}) = StringVarColHandler()
varcolhandler(::Type{<:Vector{Vector{T}}} where T) = NumericVarColHandler()
varcolhandler(::Type) = UnsupportedVarColHandler()

"""
    fits_read_var_col(f::FITSFile, colnum::Integer, data::T) where T

Read variable-length column data from a FITS file.

# Arguments
- `f`: FITS file to read from
- `colnum`: Column number (1-based)
- `data`: Pre-allocated array to store the read data
"""
function fits_read_var_col(f::FITSFile, colnum::Integer, data::T) where T
    fits_read_var_col(f, colnum, data, varcolhandler(T))
end

"""
    fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{Vector{T}}, ::NumericVarColHandler) where T

Read a variable length array column of numbers.

# Arguments
- `f`: FITS file to read from
- `colnum`: Column number (1-based)
- `data`: Pre-allocated vector of vectors to store numeric data
"""
function fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{Vector{T}}, 
                          ::NumericVarColHandler) where T
    GC.@preserve f data begin
        nrows = length(data)
        for i=1:nrows
            repeat, offset = fits_read_descript(f, colnum, i)
            data[i] = Vector{T}(undef, repeat)
            fits_read_col(f, colnum, i, 1, data[i])
        end
    end
end

"""
    fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{String}, ::StringVarColHandler)

Read a variable length array column of strings.

# Arguments
- `f`: FITS file to read from
- `colnum`: Column number (1-based)
- `data`: Pre-allocated vector to store string data
"""
function fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{String}, 
                          ::StringVarColHandler)
    status = Ref{Cint}(0)
    bufptr = Ref{Ptr{UInt8}}()
    
    GC.@preserve f data begin
        for i=1:length(data)
            repeat, offset = fits_read_descript(f, colnum, i)
            buffer = Vector{UInt8}(undef, repeat)
            
            GC.@preserve buffer begin
                bufptr[] = pointer(buffer)
                ccall((:ffgcvs, libcfitsio), Cint,
                      (Ptr{Cvoid}, Cint, Int64, Int64, Int64,
                       Ptr{UInt8}, Ref{Ptr{UInt8}}, Ptr{Cint}, Ref{Cint}),
                      f.ptr, colnum, i, 1, repeat, " ", bufptr, C_NULL, status)
                fits_assert_ok(status[])

                # Create string out of the buffer, terminating at null characters
                zeropos = something(findfirst(isequal(0x00), buffer), 0)
                data[i] = (zeropos >= 1) ? String(view(buffer, 1:(zeropos-1))) :
                                          String(buffer)
            end
        end
    end
end

"""
    fits_read_var_col(f::FITSFile, colnum::Integer, data::Any, ::UnsupportedVarColHandler)

Error handler for unsupported variable length column types.

# Arguments
- `f`: FITS file
- `colnum`: Column number
- `data`: Data of unsupported type
"""
function fits_read_var_col(f::FITSFile, colnum::Integer, data::Any, 
                          ::UnsupportedVarColHandler)
    error("Unsupported variable length column type")
end

"""
    fits_write_var_col(f::FITSFile, colnum::Integer, data::T) where T

Write variable-length column data to a FITS file.

# Arguments
- `f`: FITS file to write to
- `colnum`: Column number (1-based)
- `data`: Data to write (Vector{String} or Vector{Vector{T}} where T<:Number)
"""
function fits_write_var_col(f::FITSFile, colnum::Integer, data::T) where T
    fits_write_var_col(f, colnum, data, varcolhandler(T))
end

"""
    fits_write_var_col(f::FITSFile, colnum::Integer, data::Vector{String}, ::StringVarColHandler)

Write a variable length string column to a FITS file.

# Arguments
- `f`: FITS file to write to
- `colnum`: Column number (1-based)
- `data`: Vector of strings to write
"""
function fits_write_var_col(f::FITSFile, colnum::Integer, data::Vector{String}, 
                           ::StringVarColHandler)
    for el in data
        fits_assert_isascii(el)
    end
    
    status = Ref{Cint}(0)
    buffer = Ref{Ptr{UInt8}}()

    GC.@preserve f data begin
        for i=1:length(data)
            str_data = data[i]
            GC.@preserve str_data begin
                buffer[] = pointer(str_data)
                ccall((:ffpcls, libcfitsio), Cint,
                      (Ptr{Cvoid}, Cint, Int64, Int64, Int64, Ref{Ptr{UInt8}},
                       Ref{Cint}),
                      f.ptr, colnum, i, 1, length(str_data), buffer, status)
                fits_assert_ok(status[])
            end
        end
    end
end

"""
    fits_write_var_col(f::FITSFile, colnum::Integer, data::Vector{Vector{T}}, ::NumericVarColHandler) where T<:FITSTableScalar

Write a variable length numeric column to a FITS file.

# Arguments
- `f`: FITS file to write to
- `colnum`: Column number (1-based)
- `data`: Vector of numeric vectors to write
"""
function fits_write_var_col(f::FITSFile, colnum::Integer, data::Vector{Vector{T}}, 
                           ::NumericVarColHandler) where T<:FITSTableScalar
    GC.@preserve f data begin
        for i=1:length(data)
            row_data = data[i]
            GC.@preserve row_data begin
                fits_write_col(f, colnum, i, 1, row_data)
            end
        end
    end
end

"""
    fits_write_var_col(f::FITSFile, colnum::Integer, data::Any, ::UnsupportedVarColHandler)

Error handler for unsupported variable length column types.

# Arguments
- `f`: FITS file
- `colnum`: Column number
- `data`: Data of unsupported type

# Throws
- `ErrorException` with message about unsupported type
"""
function fits_write_var_col(f::FITSFile, colnum::Integer, data::Any, 
                           ::UnsupportedVarColHandler)
    error("column data must be a leaf type: e.g., Vector{Vector{Int}}, not Vector{Vector{T}}.")
end

"""
    fits_write_var_col(::Type{ASCIITableHDU}, ::Any)

Error handler for ASCII table HDUs, which do not support variable length columns.

# Throws
- `ErrorException` indicating lack of support for variable length columns
"""
fits_write_var_col(::Type{ASCIITableHDU}, ::Any) = 
    error("variable length columns not supported in ASCII tables")

"""
    write_internal(f::FITS, colnames::Vector{String}, coldata::Vector, hdutype, name, ver, header, units, varcols)

Add a new TableHDU to a FITS object with support for variable length columns.

# Arguments
- `f`: FITS object to write to
- `colnames`: Vector of column names
- `coldata`: Vector of column data
- `hdutype`: Type of HDU to create
- `name`: Name of the extension
- `ver`: Version number of the extension
- `header`: Optional header to write
- `units`: Optional units for columns
- `varcols`: Specification of variable length columns
"""
function write_internal(f::FITS, colnames::Vector{String},
                       coldata::Vector, hdutype, name, ver, header, units,
                       varcols)
    fits_assert_open(f.fitsfile)
    for el in colnames
        fits_assert_isascii(el)
    end

    GC.@preserve f colnames coldata begin
        # Add header if it exists
        local_header = if header isa FITSHeader
            deepcopy(header)
        else
            nothing
        end

        # move to last HDU; table will be added after the CHDU
        nhdus = Int(fits_get_num_hdus(f.fitsfile))
        (nhdus > 1) && fits_movabs_hdu(f.fitsfile, nhdus)

        ncols = length(colnames)
        ttype = [pointer(name) for name in colnames]

        # determine which columns are requested to be variable-length
        isvarcol = zeros(Bool, ncols)
        if !isnothing(varcols)
            for i=1:ncols
                isvarcol[i] = (i in varcols) || (colnames[i] in varcols)
            end
        end

        # create an array of tform strings (which we will create pointers to)
        tform_str = Vector{String}(undef, ncols)
        for i in 1:ncols
            tform_str[i] = if isvarcol[i]
                fits_tform_v(hdutype, coldata[i])
            else
                fits_tform(hdutype, coldata[i])
            end
        end
        tform = [pointer(s) for s in tform_str]

        # get units
        tunit = if isnothing(units)
            C_NULL
        else
            Ptr{UInt8}[(haskey(units, n) ? pointer(units[n]) : C_NULL)
                       for n in colnames]
        end

        # extension name
        name_ptr = isnothing(name) ? Ptr{UInt8}(C_NULL) : pointer(name)

        # Create table
        status = Ref{Cint}(0)
        ccall(("ffcrtb", libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Int64, Cint, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}},
               Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ref{Cint}),
              f.fitsfile.ptr, table_type_code(hdutype), 0, ncols,
              ttype, tform, tunit, name_ptr, status)
        fits_assert_ok(status[])

        # For binary tables, write tdim info
        if hdutype === TableHDU
            for (i, a) in enumerate(coldata)
                isvarcol[i] || fits_write_tdim(f.fitsfile, i, fits_tdim(a))
            end
        end

        # Write header if exists
        if local_header !== nothing
            fits_write_header(f.fitsfile, local_header, true)
        end
        
        if isa(ver, Integer)
            fits_update_key(f.fitsfile, "EXTVER", ver)
        end

        # Write column data
        for (i, a) in enumerate(coldata)
            if isvarcol[i]
                fits_write_var_col(f.fitsfile, i, a)
            else
                fits_write_col(f.fitsfile, i, 1, 1, a)
            end
        end
    end
    nothing
end

"""
    write(f::FITS, colnames, coldata;
          hdutype=TableHDU, name=nothing, ver=nothing,
          header=nothing, units=nothing, varcols=nothing)

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
    write(f::FITS, data::Dict;
          hdutype=TableHDU, name=nothing, ver=nothing,
          header=nothing, units=nothing, varcols=nothing)

Create a new table extension and write data to it. If the FITS file is
currently empty then a dummy primary array will be created before
appending the table extension to it. `data` should be a dictionary
with String keys (giving the column names) and Array values
(giving data to write to each column). The following types are
supported in binary tables: `UInt8`, `Int8`, `UInt16`, `Int16`,
`UInt32`, `Int32`, `UInt64`, `Int64`, `Float32`, `Float64`, `Complex{Float32}`,
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
function write(f::FITS, data::AbstractDict{<:AbstractString};
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
        data[i] = Vector{T}(undef, repeat)
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
        buffer = Vector{UInt8}(undef, repeat)
        bufptr[] = pointer(buffer)
        ccall((:ffgcvs, libcfitsio), Cint,
              (Ptr{Cvoid}, Cint, Int64, Int64, Int64,
               Ptr{UInt8}, Ref{Ptr{UInt8}}, Ptr{Cint}, Ref{Cint}),
              f.ptr, colnum, i, 1, repeat, " ", bufptr, C_NULL, status)
        fits_assert_ok(status[])

        # Create string out of the buffer, terminating at null characters
        zeropos = something(findfirst(isequal(0x00), buffer), 0)
        data[i] = (zeropos >= 1) ? String(buffer[1:(zeropos-1)]) :
                                   String(buffer)
    end
end

# Read a table column
function read(hdu::ASCIITableHDU, colname::String; case_sensitive::Bool=true)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    ### TODO: the following `if` is a deprecation warning for the case-sensitivity change.
    ### Remove it when we want to remove the deprecation.
    # `case_sensitive=true` is the default behavior, the only case where the deprecation is
    # necessary.
    if case_sensitive
        colnum = try
            fits_get_colnum(hdu.fitsfile, colname, case_sensitive=true)
        catch
            # This temporary variable is to allow showing the depwarn only if the function
            # didn't error.
            tmp = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=false)
            Base.depwarn("The new default behavior of `read(::ASCIITableHDU, ::String)`\n" *
                         "is to require case-sensitive column names.  Either pass the column name\n"*
                         "with proper capitalization or use `read(hdu, \"$(colname)\", case_sensitive=false)`",
                         :read)
            tmp
        end
    else
        colnum = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=case_sensitive)
    end

    typecode, repcnt, width = fits_get_eqcoltype(hdu.fitsfile, colnum)
    T = CFITSIO_COLTYPE[typecode]

    result = Vector{T}(undef, nrows)
    fits_read_col(hdu.fitsfile, colnum, 1, 1, result)

    return result
end


"""
    read(hdu::TableHDU, colname; case_sensitive=true)

Read a column as an array from the given table HDU.

The column name may contain wild card characters (`*`, `?`, or
`#`). The `*` wild card character matches any sequence of
characters (including zero characters) and the `?` character
matches any single character. The `#` wildcard will match any
consecutive string of decimal digits (0-9). The string must match a
unique column.  The optional boolean keyword `case_sensitive`,
`true` by default, specifies whether the column name is to be
considered case sensitive.

!!! note "Array order"

    Julia arrays are column-major (like Fortran), not row-major (like C
    and numpy), so elements of multi-dimensional columns will be the
    transpose of what you get with astropy.
"""
function read(hdu::TableHDU, colname::String; case_sensitive::Bool=true)
    assert_open(hdu)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    ### TODO: the following `if` is a deprecation warning for the case-sensitivity change.
    ### Remove it when we want to remove the deprecation.
    # `case_sensitive=true` is the default behavior, the only case where the deprecation is
    # necessary.
    if case_sensitive
        colnum = try
            fits_get_colnum(hdu.fitsfile, colname, case_sensitive=true)
        catch
            # This temporary variable is to allow showing the depwarn only if the function
            # didn't error.
            tmp = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=false)
            Base.depwarn("The new default behavior of `read(::TableHDU, ::String)`\n" *
                         "is to require case-sensitive column names.  Either pass the column name\n"*
                         "with proper capitalization or use `read(hdu, \"$(colname)\", case_sensitive=false)`",
                         :read)
            tmp
        end
    else
        colnum = fits_get_colnum(hdu.fitsfile, colname, case_sensitive=case_sensitive)
    end


    T, rowsize, isvariable = fits_get_col_info(hdu.fitsfile, colnum)

    result = Array{T}(undef, rowsize..., nrows)

    if isvariable
        fits_read_var_col(hdu.fitsfile, colnum, result)
    else
        fits_read_col(hdu.fitsfile, colnum, 1, 1, result)
    end

    return result
end

#Tables.jl integration

const EitherTableHDU = Union{TableHDU, ASCIITableHDU}
Tables.istable(::Type{<:EitherTableHDU}) = true
Tables.columnaccess(::Type{<:EitherTableHDU}) = true
Tables.columns(t::EitherTableHDU) = t

function Tables.columnnames(t::EitherTableHDU)
    cns = FITSIO.colnames(t)
    #filter out variable-length cols
    isvar = [last(fits_get_col_info(t.fitsfile, i)) for i in 1:length(cns)]
    Symbol.(cns[.! isvar])
end
function Tables.getcolumn(t::EitherTableHDU, s::Symbol)
    col = FITSIO.read(t, String(s))
    if fits_get_col_info(t.fitsfile, findfirst(FITSIO.colnames(t) .== String(s)))[3]
        error("variable-length columns not supported in Tables.jl interface")
    end
    #reshape multidimensional array into array of (possibly multidimensional) arrays if necessary
    if (dim = length(size(col))) == 1
        col
    else
        [col[vcat(repeat([:], dim-1), i)...] for i in 1:size(col, dim)]
    end
end
#this uses Tables.columnnames rather that FITSIO.colnames so as to ignore variable length columns
Tables.getcolumn(t::EitherTableHDU, i::Int) = Tables.getcolumn(t, Tables.columnnames(t)[i])

