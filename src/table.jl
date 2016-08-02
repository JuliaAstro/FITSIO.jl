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
                         (Compat.ASCIIString, 'A',  16),
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
typealias FITSTableScalar @compat(Union{UInt8, Int8, Bool, UInt16, Int16, UInt32,
                                Int32, Int64, Float32, Float64, Complex64,
                                Complex128})

# Helper function for reading information about a (binary) table column
# Returns: (eltype, rowsize, isvariable)
function fits_get_col_info(f::FITSFile, colnum::Integer)
    eqtypecode, repeat, width = fits_get_eqcoltype(f, colnum)
    isvariable = eqtypecode < 0
    eqtypecode = abs(eqtypecode)

    (eqtypecode == 1) && error("BitArray ('X') columns not yet supported")

    T = CFITSIO_COLTYPE[eqtypecode]

    if isvariable
        if T !== Compat.ASCIIString
            T = Vector{T}
        end
        rowsize = Int[]
    else
        if T === Compat.ASCIIString
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
function var_col_maxlen(tform::Compat.ASCIIString)
    maxlen = -1
    i = search(tform, '(')
    if i > 0
        j = search(tform, ')', i)
        if j > 0
            try maxlen = parseint(tform[i+1:j-1]) end
        end
    end
    return maxlen
end

# Helper function for getting fits tdim shape for given array
fits_tdim(A::Array) = (ndims(A) == 1)? [1]: [size(A, i) for i=1:ndims(A)-1]
function fits_tdim(A::Array{Compat.ASCIIString})
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
fits_tform(::Type{TableHDU}, A::Vector{Compat.ASCIIString}) = "$(maximum(length, A))A"
fits_tform(::Type{TableHDU}, A::Array{Compat.ASCIIString}) = "$(prod(fits_tdim(A)))A$(maximum(length, A))"

# variable length columns
fits_tform_v{T<:FITSTableScalar}(::Type{TableHDU}, A::Vector{Vector{T}}) = "1P$(fits_tform_char(T))($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{Compat.ASCIIString}) = "1PA($(maximum(length(A))))"
fits_tform_v(::Type{TableHDU}, A::Vector{Vector}) = error("column data must be a leaf type: e.g., Vector{Vector{Int}}, not Vector{Vector{T}}.")
fits_tform_v(::Type{TableHDU}, ::Any) = error("variable length columns only supported for arrays of arrays and arrays of Compat.ASCIIString")
fits_tform_v(::Type{ASCIITableHDU}, ::Any) = error("variable length columns not supported in ASCII tables")


fits_tform(::Type{ASCIITableHDU}, ::Vector{Int16}) = "I7"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Int32}) = "I12"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float32}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, ::Vector{Float64}) = "E26.17"
fits_tform(::Type{ASCIITableHDU}, A::Vector{Compat.ASCIIString}) = "A$(maximum(length, A))"
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
    colnames_in = [Array(UInt8, 70) for i=1:ncols]
    coltforms_in = [Array(UInt8, 70) for i=1:ncols]
    nrows = Array(Int64, 1)
    status = Cint[0]

    # fits_read_btblhdrll (Can pass NULL for return fields not needed.)
    ccall(("ffghbnll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Cint}, Ptr{Ptr{UInt8}},  # nrows, tfields, ttype
           Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}}, Ptr{UInt8},  # tform,tunit,extname
           Ptr{Clong}, Ptr{Cint}),  # pcount, status
          hdu.fitsfile.ptr, ncols, nrows, C_NULL, colnames_in, coltforms_in,
          C_NULL, C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    colnames = [unsafe_string(pointer(item)) for item in colnames_in]
    coltforms = [unsafe_string(pointer(item)) for item in coltforms_in]

    # get some more information for all the columns
    coltypes = Array(Compat.ASCIIString, ncols)
    colrowsizes = Array(Compat.ASCIIString, ncols)
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
    Rows: $(nrows[1])
    Columns: """)
    show_ascii_table(
        io, ["Name", "Size", "Type", "TFORM"],
        Vector{Compat.ASCIIString}[colnames, colrowsizes, coltypes, coltforms],
        2, 9)
    if showlegend
        print(io, "\n         (*) = variable-length column\n")
    end
end

function show(io::IO, hdu::ASCIITableHDU)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)
    ncols = fits_get_num_cols(hdu.fitsfile)

    # allocate return arrays for column names & types
    colnames_in = [Array(UInt8, 70) for i=1:ncols]
    coltforms_in = [Array(UInt8, 70) for i=1:ncols]
    nrows = Array(Int64, 1)
    status = Cint[0]

    # fits_read_atblhdrll (Can pass NULL for return fields not needed)
    ccall(("ffghtbll", libcfitsio), Cint,
          (Ptr{Void}, Cint,  # Inputs: fitsfile, maxdim
           Ptr{Int64}, Ptr{Int64}, Ptr{Cint},  # rowlen, nrows, tfields
           Ptr{Ptr{UInt8}}, Ptr{Clong}, Ptr{Ptr{UInt8}},  # ttype, tbcol, tform
           Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ptr{Cint}),  # tunit, extname, status
          hdu.fitsfile.ptr, ncols,
          C_NULL, nrows, C_NULL,
          colnames_in, C_NULL, coltforms_in,
          C_NULL, C_NULL, status)
    fits_assert_ok(status[1])

    # parse out results
    colnames = [unsafe_string(pointer(item)) for item in colnames_in]
    coltforms = [unsafe_string(pointer(item)) for item in coltforms_in]

    # Get additional info
    coltypes = Array(Compat.ASCIIString, ncols)
    for i in 1:ncols
        eqtypecode, repeat, width = fits_get_eqcoltype(hdu.fitsfile, i)
        T = CFITSIO_COLTYPE[eqtypecode]
        coltypes[i] = repr(T)        
    end

    print(io, """
    File: $(fits_file_name(hdu.fitsfile))
    HDU: $(hdu.ext)$(fits_get_ext_info_string(hdu.fitsfile))
    Type: ASCIITable
    Rows: $(nrows[1])
    Columns: """)
    show_ascii_table(io, ["Name", "Type", "TFORM"],
                     Vector{Compat.ASCIIString}[colnames, coltypes, coltforms], 2, 9)
end

# Write a variable length array column of numbers
# (separate implementation from normal fits_write_col function because
#  we must make separate calls to `fits_write_col` for each row.)
function fits_write_var_col{T}(f::FITSFile, colnum::Integer,
                               data::Vector{Vector{T}})
    for i=1:length(data)
        fits_write_col(f, colnum, i, 1, data[i])
    end
end

function fits_write_var_col(f::FITSFile, colnum::Integer,
                            data::Vector{Compat.ASCIIString})
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
              (Ptr{Void}, Cint, Int64, Int64, Int64, Ptr{Ptr{UInt8}},
               Ptr{Cint}),
              f.ptr, colnum, i, 1, length(data[i]), buffer, status)
        fits_assert_ok(status[1])
    end
end

# Add a new TableHDU to a FITS object
function write_internal(f::FITS, colnames::Vector{Compat.ASCIIString},
                        coldata::Vector, hdutype, name, ver, header, units,
                        varcols)
    fits_assert_open(f.fitsfile)

    # move to last HDU; table will be added after the CHDU
    nhdus = @compat(Int(fits_get_num_hdus(f.fitsfile)))
    (nhdus > 1) && fits_movabs_hdu(f.fitsfile, nhdus)

    ncols = length(colnames)
    ttype = [pointer(name) for name in colnames]

    # determine which columns are requested to be variable-length
    isvarcol = zeros(Bool, ncols)
    if !isa(varcols, @compat(Void))
        for i=1:ncols
            isvarcol[i] = (i in varcols) || (colnames[i] in varcols)
        end
    end

    # create an array of tform strings (which we will create pointers to)
    tform_str = Array(Compat.ASCIIString, ncols)
    for i in 1:ncols
        if isvarcol[i]
            tform_str[i] = fits_tform_v(hdutype, coldata[i])
        else
            tform_str[i] = fits_tform(hdutype, coldata[i])
        end
    end
    tform = [pointer(s) for s in tform_str]

    # get units
    if isa(units, @compat(Void))
        tunit = C_NULL
    else
        tunit = Ptr{UInt8}[(haskey(units, n)? pointer(units[n]): C_NULL)
                           for n in colnames]
    end

    # extension name
    name_ptr = (isa(name, @compat(Void)) ? convert(Ptr{UInt8}, C_NULL) :
                   pointer(name))

    status = Cint[0]
    ccall(("ffcrtb", libcfitsio), Cint,
          (Ptr{Void}, Cint, Int64, Cint, Ptr{Ptr{UInt8}}, Ptr{Ptr{UInt8}},
           Ptr{Ptr{UInt8}}, Ptr{UInt8}, Ptr{Cint}),
          f.fitsfile.ptr, table_type_code(hdutype), 0, ncols,  # 0 = nrows
          ttype, tform, tunit, name_ptr, status)
    fits_assert_ok(status[1])

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

function write(f::FITS, colnames::Vector{Compat.ASCIIString}, coldata::Vector;
               units=nothing, header=nothing, hdutype=TableHDU,
               name=nothing, ver=nothing, varcols=nothing)
    if length(colnames) != length(coldata)
        error("length of colnames and length of coldata must match")
    end
    write_internal(f, colnames, coldata, hdutype, name, ver, header,
                   units, varcols)
end

function write(f::FITS, data::Dict{Compat.ASCIIString};
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
function fits_read_var_col{T}(f::FITSFile, colnum::Integer,
                              data::Vector{Vector{T}})
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
function fits_read_var_col(f::FITSFile, colnum::Integer, data::Vector{Compat.ASCIIString})
    status = Cint[0]
    bufptr = Array(Ptr{UInt8}, 1)  # holds a pointer to the current row buffer
    for i=1:length(data)
        repeat, offset = fits_read_descript(f, colnum, i)
        buffer = Array(UInt8, repeat)
        bufptr[1] = pointer(buffer)
        ccall((:ffgcvs, libcfitsio), Cint,
              (Ptr{Void}, Cint, Int64, Int64, Int64,
               Ptr{UInt8}, Ptr{Ptr{UInt8}}, Ptr{Cint}, Ptr{Cint}),
              f.ptr, colnum, i, 1, repeat, " ", bufptr, C_NULL, status)
        fits_assert_ok(status[1])

        # Create string out of the buffer, terminating at null characters
        zeropos = search(buffer, 0x00)
        data[i] = (zeropos >= 1) ? Compat.ASCIIString(buffer[1:(zeropos-1)]) :
                                   Compat.ASCIIString(buffer)
    end
end

# Read a table column
function read(hdu::ASCIITableHDU, colname::Compat.ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname)

    typecode, repcnt, width = fits_get_eqcoltype(hdu.fitsfile, colnum)
    T = CFITSIO_COLTYPE[typecode]

    result = Array(T, nrows)
    fits_read_col(hdu.fitsfile, colnum, 1, 1, result)

    return result
end

function read(hdu::TableHDU, colname::Compat.ASCIIString)
    fits_assert_open(hdu.fitsfile)
    fits_movabs_hdu(hdu.fitsfile, hdu.ext)

    nrows = fits_get_num_rows(hdu.fitsfile)
    colnum = fits_get_colnum(hdu.fitsfile, colname)

    T, rowsize, isvariable = fits_get_col_info(hdu.fitsfile, colnum)

    result = Array(T, rowsize..., nrows)

    if isvariable
        fits_read_var_col(hdu.fitsfile, colnum, result)
    else
        fits_read_col(hdu.fitsfile, colnum, 1, 1, result)
    end

    return result
end
